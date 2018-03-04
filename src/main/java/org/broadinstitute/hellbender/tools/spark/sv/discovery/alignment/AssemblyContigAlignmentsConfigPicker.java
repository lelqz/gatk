package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoverFromLocalAssemblyContigAlignmentsSpark;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.SvDiscoveryUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVUtils;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;

public class AssemblyContigAlignmentsConfigPicker {

    /**
     * A filter that is used to remove contigs upfront which doesn't meet the following criteria
     * either:
     *  has only 1 mapping, with MQ strictly above this threshold
     * or:
     *  has more than 1 mappings, but only 1 mapping has MQ strictly above this threshold and it has a large gap in it.
     */
    static final int ALIGNMENT_MQ_THRESHOLD = 20;

    /**
     * A filter to boost configuration scoring implemented here:
     * if the configuration has more than 10 mappings, then
     * any mappings in such configuration with MQ
     * not strictly above this threshold is classified as bad and filtered.
     */
    static final int ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST = 10;

    /**
     * Filters an input of SAM file containing alignments of a single-ended long read that
     * aims at providing an "optimal coverage" of the assembly contig, based on an heuristic scoring scheme
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     *
     * @param assemblyAlignments    long read alignments
     * @param header                header for the long reads
     * @param scoreDiffTolerance    a tolerance where if two configurations' scores differ by less than or equal to this amount, they are considered equally good
     * @param toolLogger            logger for, most likely, debugging uses
     *
     * @return              contigs with alignments filtered and custom formatted as {@link AlignmentInterval}
     */
    // TODO: 3/1/18 change interface here: not to return AlignedContig, but return AssemblyContigWithFineTunedAlignments
    public static JavaRDD<AlignedContig> createOptimalCoverageAlignmentSetsForContigs(final JavaRDD<GATKRead> assemblyAlignments,
                                                                                      final SAMFileHeader header,
                                                                                      final String nonCanonicalContigNamesFile,
                                                                                      final Double scoreDiffTolerance,
                                                                                      final Logger toolLogger) {

        final JavaRDD<AlignedContig> parsedContigAlignments =
                convertRawAlignmentsToAlignedContigAndFilterByQuality(assemblyAlignments, header, toolLogger);

        return filterAndSplitGappedAlignmentInterval(parsedContigAlignments, nonCanonicalContigNamesFile,
                                      header.getSequenceDictionary(), scoreDiffTolerance);
    }

    //==================================================================================================================

    /**
     * Parses input alignments into custom {@link AlignmentInterval} format, and
     * performs a primitive filtering implemented in
     * {@link #notDiscardForBadMQ(AlignedContig)} that
     * gets rid of contigs with no good alignments.
     */
    // TODO: 3/1/18 change interface here: not to return AlignedContig, but return AssemblyContigWithFineTunedAlignments
    private static JavaRDD<AlignedContig> convertRawAlignmentsToAlignedContigAndFilterByQuality(final JavaRDD<GATKRead> assemblyAlignments,
                                                                                                final SAMFileHeader header,
                                                                                                final Logger toolLogger) {
        assemblyAlignments.cache();
        toolLogger.info( "Processing " + assemblyAlignments.count() + " raw alignments from " +
                         assemblyAlignments.map(GATKRead::getName).distinct().count() + " contigs.");

        final JavaRDD<AlignedContig> parsedContigAlignments =
                new SvDiscoverFromLocalAssemblyContigAlignmentsSpark.SAMFormattedContigAlignmentParser(assemblyAlignments, header, false)
                        .getAlignedContigs()
                        .filter(AssemblyContigAlignmentsConfigPicker::notDiscardForBadMQ).cache();
        assemblyAlignments.unpersist();
        toolLogger.info( "Filtering on MQ left " + parsedContigAlignments.count() + " contigs.");
        return parsedContigAlignments;
    }

    /**
     * Idea is to keep mapped contig that
     *  either has at least two alignments over {@link #ALIGNMENT_MQ_THRESHOLD},
     *  or in the case of a single alignment, it must be MQ > {@link #ALIGNMENT_MQ_THRESHOLD}.
     * Note that we are not simply filtering out contigs with only 1 alignment because
     * they might contain large (> 50) gaps hence should be kept.
     *
     * a point that could use improvements:
     *   the current implementation exhaustively checks the power set of all possible alignments of each assembly contig,
     *   which is computationally impossible for contigs having many-but-barely-any-good alignments, yet bringing in no value,
     *   hence this primitive filtering step to get rid of these bad assembly contigs.
     */
    private static boolean notDiscardForBadMQ(final AlignedContig contig) {
        if (contig.alignmentIntervals.size() < 2 ) {
            return (!contig.alignmentIntervals.isEmpty()) && contig.alignmentIntervals.get(0).mapQual > ALIGNMENT_MQ_THRESHOLD;
        } else {
            // TODO: 3/1/18 a bug is present here that even though only one alignment has not-bad MQ, it could contain a large gap, currently it is being filtered away;
            //      we should keep the single not-bad mapping and mark the others as bad;
            //      note that a follow up fix in AssemblyContigAlignmentSignatureClassifier to classify such contigs not as incomplete
            return contig.alignmentIntervals.stream().mapToInt(ai -> ai.mapQual).filter(mq -> mq > ALIGNMENT_MQ_THRESHOLD).count() > 1;
        }
    }

    //==================================================================================================================

    /**
     * For each assembly contig, scores its alignment configurations and pick the best one(s),
     * then reconstruct the contig's alignment configuration through {@link #reConstructContigFromPickedConfiguration(Tuple2)}.
     *
     * Note that this step is essentially a flatMap operation, meaning that one contig may return 1 or multiple contigs:
     *  *) when 1 contig is yielded, it means the contig has only 1 configuration that scored the best
     *  *) when multiple contigs are yielded, it means the contig has several top-scored configurations
     * How to handle the second scenario can be treated in a separate logic unit.
     */
    @VisibleForTesting
    static JavaRDD<AlignedContig> filterAndSplitGappedAlignmentInterval(final JavaRDD<AlignedContig> parsedContigAlignments,
                                                                        final String nonCanonicalContigNamesFile,
                                                                        final SAMSequenceDictionary dictionary,
                                                                        final Double scoreDiffTolerance) {

        final Set<String> canonicalChromosomes = SvDiscoveryUtils.getCanonicalChromosomes(nonCanonicalContigNamesFile, dictionary);

        return parsedContigAlignments
                .mapToPair(alignedContig -> new Tuple2<>(new Tuple2<>(alignedContig.contigName,alignedContig.contigSequence),
                                pickBestConfigurations(alignedContig, canonicalChromosomes, scoreDiffTolerance)))
                .map(nameSeqAndConfigurations -> new Tuple2<>(nameSeqAndConfigurations._1, breakTieByPreferringLessAlignments(nameSeqAndConfigurations._2, 0)))
                .map(nameSeqAndConfigurations -> new Tuple2<>(nameSeqAndConfigurations._1, nameSeqAndConfigurations._2.stream().map(r -> r.goodMappings).collect(Collectors.toList()))) // todo temporary, later need to change interface of the next map
                .flatMap(AssemblyContigAlignmentsConfigPicker::reConstructContigFromPickedConfiguration);
    }

    /**
     * After configuration scoring and picking, the original alignments can be classified as
     * good and bad mappings:
     * good: the ones present the picked configuration
     * bad: the brings more noise than information,
     *      they can be turned into string representation following the format as in {@link AlignmentInterval#toPackedString()}
     */
    static final class GoodAndBadMappings {

        final List<AlignmentInterval> goodMappings;
        final List<AlignmentInterval> badMappings;

        GoodAndBadMappings(final List<AlignmentInterval> goodMappings, final List<AlignmentInterval> badMappings) {
            this.goodMappings = goodMappings;
            this.badMappings = badMappings;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            final GoodAndBadMappings that = (GoodAndBadMappings) o;

            if (!goodMappings.equals(that.goodMappings)) return false;
            return badMappings.equals(that.badMappings);
        }

        @Override
        public int hashCode() {
            int result = goodMappings.hashCode();
            result = 31 * result + badMappings.hashCode();
            return result;
        }
    }

    /**
     * Pick the best configurations based on a heuristic scoring scheme implemented in
     * {@link #computeScoreOfConfiguration(List, Set, int)}.
     * @return a 2-D list, where in the case when multiple configurations are equally top-scored, all such configurations are picked up
     */
    @VisibleForTesting
    static List<GoodAndBadMappings> pickBestConfigurations(final AlignedContig alignedContig,
                                                           final Set<String> canonicalChromosomes,
                                                           final Double scoreDiffTolerance) {

        // group 1: get max aligner score of mappings to canonical chromosomes and speed up in case of too many mappings
        final int maxCanonicalChrAlignerScore = alignedContig.alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // speed up if number of alignments is too high (>10)
        // if mapped to canonical chromosomes, MQ must be >10; otherwise, must have AS higher than max canonical aligner score
        // TODO: 3/1/18 should keep the bad ones and save them in the returned value
        final List<AlignmentInterval> alignmentIntervals;
        if (alignedContig.alignmentIntervals.size() > 10) {
            alignmentIntervals = alignedContig.alignmentIntervals.stream()
                    .filter(alignmentInterval -> (!canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())
                                                                           && alignmentInterval.alnScore > maxCanonicalChrAlignerScore)
                                                 || alignmentInterval.mapQual > ALIGNMENT_MQ_THRESHOLD_FOR_SPEED_BOOST)
                    .collect(Collectors.toList());
        } else {
            alignmentIntervals = alignedContig.alignmentIntervals;
        }

        final int newMaxCanonicalChrAlignerScore = alignmentIntervals.stream()
                .filter(alignmentInterval -> canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig()))
                .mapToInt(ai -> ai.alnScore).max().orElse(0); // possible that no mapping to canonical chromosomes

        // group 2: generate, and score configurations
        final List<List<AlignmentInterval>> allConfigurations = Sets.powerSet(new HashSet<>(alignmentIntervals))
                .stream().map(ArrayList::new)
                // make sure within each configuration, alignments would be sorted as they would be in a corresponding AlignedContig
                .map(ls -> ls.stream().sorted(AlignedContig.getAlignmentIntervalComparator()).collect(Collectors.toList()))
                .collect(Collectors.toList());

        final List<Double> scores = allConfigurations.stream()
                .map(configuration -> computeScoreOfConfiguration(configuration, canonicalChromosomes, newMaxCanonicalChrAlignerScore))
                .collect(SVUtils.arrayListCollector(allConfigurations.size()));

        // group 3: pick the best-scored configuration(s) (if multiple configurations have equally good scores, return all of them)
        final double maxScore = scores.stream().mapToDouble(Double::doubleValue).max()
                .orElseThrow(() -> new GATKException("Cannot find best-scoring configuration on alignments of contig: " + alignedContig.contigName));

        return IntStream.range(0, allConfigurations.size())
                .filter(i -> {
                    final Double s = scores.get(i);
                    // two configurations with would-be-same-scores can differ by a tolerance due to the sin of comparing
                    // could-be-close floating point values
                    // (see http://www.cygnus-software.com/papers/comparingfloats/Comparing%20floating%20point%20numbers.htm)
                    final Double tol = Math.max(Math.ulp(s), scoreDiffTolerance);
                    return s >= maxScore || maxScore - s <= tol;
                })
                .mapToObj(p -> {
                    final ArrayList<AlignmentInterval> copy = new ArrayList<>(alignmentIntervals);
                    final List<AlignmentInterval> pickedAlignments = allConfigurations.get(p);
                    copy.removeAll(pickedAlignments); // remove picked, left are bad
                    return new GoodAndBadMappings(pickedAlignments, copy);
                })
                .collect(Collectors.toList());
    }

    /**
     * quick and dirty implementation for computing score of given configuration of alignments,
     * no assumption on the ordering of input alignments
     */
    @VisibleForTesting
    static double computeScoreOfConfiguration(final List<AlignmentInterval> configuration,
                                              final Set<String> canonicalChromosomes,
                                              final int maxCanonicalChrAlignerScore) {

        final double tigExplainQual = computeTigExplainQualOfOneConfiguration(configuration, canonicalChromosomes, maxCanonicalChrAlignerScore);

        int redundancy = 0;
        for (int i = 0; i < configuration.size() -1 ; ++i) {
            for (int j = i + 1; j < configuration.size(); ++j) {
                final int overlap = AlignmentInterval.overlapOnContig(configuration.get(i), configuration.get(j));
                redundancy += overlap;
            }
        }

        return tigExplainQual - redundancy;
    }

    private static double computeTigExplainQualOfOneConfiguration(final List<AlignmentInterval> configuration,
                                                                  final Set<String> canonicalChromosomes,
                                                                  final int maxCanonicalChrAlignerScore) {
        double tigExplainedQual = 0;
        for (final AlignmentInterval alignmentInterval : configuration) {
            final int len = alignmentInterval.endInAssembledContig - alignmentInterval.startInAssembledContig + 1;
            final double weight;
            if (canonicalChromosomes.contains(alignmentInterval.referenceSpan.getContig())) {
                weight = alignmentInterval.mapQual/60.0;
            } else {
                weight = Math.max(alignmentInterval.mapQual/60.0,
                                  alignmentInterval.alnScore > maxCanonicalChrAlignerScore ? 1 : 0);
            }
            tigExplainedQual += weight * len;
        }
        return tigExplainedQual;
    }

    //==================================================================================================================

    /**
     * Reconstructs (possibly more than one) {@link AlignedContig} based on
     * the given best-scored configuration(s) in {@code nameSeqAndBestConfigurationsOfOneRead}.
     *
     * @param nameSeqAndBestConfigurationsOfOneRead the name, sequence, and picked best alignment configuration(s) of an assembly contig
     * @return The number of returned contigs will be the same as the given best-scored configurations.
     */
    private static Iterator<AlignedContig> reConstructContigFromPickedConfiguration(
            final Tuple2<Tuple2<String, byte[]>, List<List<AlignmentInterval>>> nameSeqAndBestConfigurationsOfOneRead) {

        final String contigName = nameSeqAndBestConfigurationsOfOneRead._1._1;
        final byte[] contigSeq = nameSeqAndBestConfigurationsOfOneRead._1._2;
        final List<List<AlignmentInterval>> bestConfigurations = nameSeqAndBestConfigurationsOfOneRead._2;
        if (bestConfigurations.size() > 1) { // more than one best configuration
            return bestConfigurations.stream()
                    .map(configuration ->
                            new AlignedContig(contigName, contigSeq, splitGaps(configuration),
                                    true))
                    .sorted(sortConfigurations())
                    .iterator();
        } else {
            return Collections.singletonList(
                    new AlignedContig(contigName, contigSeq, splitGaps(bestConfigurations.get(0)),
                            false))
                    .iterator();
        }
    }

    /**
     * when two configurations are the same,
     * put the one with less alignments,
     * or less summed mismatches if still tie
     * first
     */
    private static Comparator<AlignedContig> sortConfigurations() {
        final Comparator<AlignedContig> numFirst
                = (AlignedContig x, AlignedContig y) -> Integer.compare(x.alignmentIntervals.size(), y.alignmentIntervals.size());
        final Comparator<AlignedContig> mismatchSecond
                = (AlignedContig x, AlignedContig y) -> Integer.compare(x.alignmentIntervals.stream().mapToInt(ai -> ai.mismatches).sum(),
                                                                        y.alignmentIntervals.stream().mapToInt(ai -> ai.mismatches).sum());
        return numFirst.thenComparing(mismatchSecond);
    }

    @VisibleForTesting
    static List<AlignmentInterval> splitGaps(final List<AlignmentInterval> configuration) {

        // TODO: 3/1/18 don't ditch the bad ones, instead save them, but mark accordingly
        // 1st pass, split gapped alignments when available
        final List<Iterable<AlignmentInterval>> alignmentSplitChildren =
                configuration.stream()
                        .map(alignment -> {
                            final Iterable<AlignmentInterval> split;
                            if (alignmentContainsLargeGap(alignment)) {
                                split = ContigAlignmentsModifier.splitGappedAlignment(alignment, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                                        SvCigarUtils.getUnclippedReadLength(alignment.cigarAlong5to3DirectionOfContig));
                            } else {
                                split = Collections.singletonList(alignment);
                            }
                            return split;
                        }).collect(Collectors.toList());

        // 2nd pass make a choice between gapped and overlapping alignment (alignments that are not favored has its "2nd" set to null)
        final int count = configuration.size();
        for (int i = 0; i < count; ++i) {
            final AlignmentInterval alignment = configuration.get(i);
            final Iterable<AlignmentInterval> split = alignmentSplitChildren.get(i);
            if ( split != null && Iterables.size(split) != 1 ) { // the split could be null (i.e. to be filtered out), or could contain no large gaps (i.e. should not check it)
                for (int j = 0; j < count; ++j) {
                    final AlignmentInterval other = configuration.get(j);
                    if (j == i || AlignmentInterval.overlapOnContig(alignment, other) == 0)
                        continue;

                    if ( Utils.stream(split).anyMatch(other::containsOnRead) ) {
                        if ( gappedAlignmentOffersBetterCoverage(alignment, other) ) {
                            alignmentSplitChildren.set(j, null);
                        } else {
                            alignmentSplitChildren.set(i, null);
                        }
                    }
                }
            }
        }

        // filter out to-be-removed and done
        return alignmentSplitChildren.stream()
                .filter(Objects::nonNull)
                .flatMap(Utils::stream).collect(Collectors.toList());
    }

    private static boolean alignmentContainsLargeGap(final AlignmentInterval alignment) {
        return alignment.cigarAlong5to3DirectionOfContig.getCigarElements().stream()
                .anyMatch(cigarElement ->
                        cigarElement.getOperator().isIndel() && cigarElement.getLength() >= GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY);
    }

    private static boolean gappedAlignmentOffersBetterCoverage(final AlignmentInterval gapped,
                                                               final AlignmentInterval overlappingNonGapped) {
        final int diff = gapped.getSizeOnRead() - overlappingNonGapped.getSizeOnRead();
        if ( diff == 0) {
            return gapped.alnScore > overlappingNonGapped.alnScore;
        } else {
            return diff > 0;
        }
    }

    //==================================================================================================================

    /**
     * For contigs with more than 1 best-scored configurations as determined by
     * {@link #pickBestConfigurations(AlignedContig, Set, Double)},
     * save the contigs that has one and only one configuration that
     * has all mapping quality strictly above the specified {@code threshold}.
     * Example:
     *  if a contig has two equal scored configurations with MQ's {10, 60, 60}, and {60, 60},
     *  this function will favor/pick the {60, 60} configuration if the threshold is 10,
     *  hence remove the ambiguity.
     */
    @VisibleForTesting
    static List<GoodAndBadMappings> breakTieByPreferringLessAlignments(
            final List<GoodAndBadMappings> differentRepresentationsForOneContig,
            final int threshold) {
        if ( differentRepresentationsForOneContig.size() == 1) {
            return differentRepresentationsForOneContig;
        } else {
            return
                    Utils.stream(differentRepresentationsForOneContig)
                            .filter( rep -> rep.goodMappings.stream().mapToInt(ai -> ai.mapQual).min().orElse(threshold) > threshold )
                            .collect(Collectors.toList());
        }
    }

    // TODO: 3/1/18 implement another simple fix that classify contigs with the following signature as ambiguous
    // one configuration with single good mapping to non-canonical chromosome,
    // another configuration with split alignments to canonical chromosomes
}
