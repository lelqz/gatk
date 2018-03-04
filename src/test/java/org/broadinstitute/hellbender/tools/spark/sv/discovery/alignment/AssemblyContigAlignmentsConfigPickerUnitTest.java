package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import com.google.common.collect.Lists;
import htsjdk.samtools.TextCigarCodec;
import org.apache.commons.collections4.IterableUtils;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.CigarUtils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SVTestUtils.fromPrimarySAMRecordString;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SVTestUtils.makeDummySequence;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.b38_canonicalChromosomes;
import static org.broadinstitute.hellbender.tools.spark.sv.discovery.SimpleSVDiscoveryTestDataProvider.b38_seqDict;
import static org.testng.Assert.assertEquals;
import static org.testng.Assert.assertTrue;

public class AssemblyContigAlignmentsConfigPickerUnitTest extends GATKBaseTest {

    @DataProvider(name = "contigAlignmentsHeuristicFilter")
    private Object[][] createTestData() {

        final List<Object[]> data = new ArrayList<>(20);

        AlignmentInterval intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 1948156, 1948936),
                1, 787, TextCigarCodec.decode("257M4I182M2I342M361S"), true, 60, 8, 733, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval intervalTwo = new AlignmentInterval(new SimpleInterval("chr21", 1948935, 1949190),
                893, 1148, TextCigarCodec.decode("892H256M"), true, 60, 3, 241, ContigAlignmentsModifier.AlnModType.NONE);
        AlignedContig contig = new AlignedContig("asm000063:tig00003", "CCACTGTGCCCGGCCAAGGGTCCCCGGTTCTGAAAGTGGAAGGGGTGCGGCTGCCTCAGGAGTCACCACGGCAACAAGAACCTGGACCTGAGCGCAGGTGGTCAGATTCTGGGGCCAGCAGCTTTTTGGTTTTTAGAGACGAGGTCTCACTCTGTTGCCCAGGCTGGAGTGCAGTGGTGCGATCACTGCACCCTGCAGCCTCGGCCTCCTGGTTTCAAGTGACCACAGATGCATGCAGCCATGCTTGGCATATATAAATATATATATATATATATATTTATGTGTATATTGGTAGAGACATGGTCTTGTTATATTGCCCAGGCTGATCGCAAACATCTGCTTAAGCGATCCTCCTGCGTTGGCTCTCCAAAGTATTGGGATTATAGGCATGAGCTACCATGGCCTGGCCTCCTTATTCTAGTCTTTTCTTTCCTTTCTTCTTGTTTTTTTTTTTTTTTTGGCAGGGTCTCACTCTGTCACCCAGGCTGCAGTGCAGTGGTGTGATCACAGCTCACTGCAGCCTCAACTTCCCAGGCTCAAGCGATCCTCCCGGCTCAGCATCCTGAGTAGCTGGGACTACAGATGCATGTCACCACGCCTGGCTAAATTTTCTTCTTTGTAGATATGGGGTCTCACCATGTAGTACTTTTCAATGTATTAAGCATCCTTATTTGATATTTGATGCCTGATAATACCCATGTCTGAACCATGCAAGATTGCTGCAATTCCTTCCTTCCTTCCCTCCCTCCTTCCCTTCCTTCCTTCCCTTTCCTTCCTTCCTCTTTCCCTCCCTTCTTTCCTTCCCTTTCCCTCCCTCCCTTCCTTCCTCTTTCCTTCCTTCCTTTCCCTCCCTTACTCCTTCCTTCCCTTCCCCTTCCTTCTTCCTTCTCTCCCTCCCTCCCTTCCCCTCCCTTACTCCCTTCCTTCCTCCTTCCCTCCCTCCTTTCCTTCATTCCCTTCCTTCCCCTTCCCCTTCCTTCCTTCTCTCCCTCCCTCCTTCCTTCCCTCCTTTCCTTCCTTCCTTCCTTTCCTTTCCCTCCTTCCTCCCTCCCTCCTTTCCTTCCTTCCTTTCCTTTCCTCCCTTCCCTCCCTCCCTCCCTCCCTTCCTTCCCCTCCCTCCCTCCTTTCCTTCTTTCGACAGAGTCTTG".getBytes(),
                Arrays.asList(intervalOne, intervalTwo), false);
        data.add(new Object[]{contig, Arrays.asList(intervalOne), Arrays.asList(intervalOne, intervalTwo), 1, 2});

        intervalOne = new AlignmentInterval(new SimpleInterval("chr2", 1422222, 1422435),
                1, 270,  TextCigarCodec.decode("75M56I139M"), false, 60, 56, 142, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr2_KI270774v1_alt", 105288, 105557),
                1, 270,  TextCigarCodec.decode("114M1I27M1I127M"), false, 56, 13, 179, ContigAlignmentsModifier.AlnModType.NONE);
        contig = new AlignedContig("asm002608:tig00001", "ATGCTGGGGAATTTGTGTGCTCCTTGGGTGGGGACGAGCATGGAAGGCGCGTGGGACTGAAGCCTTGAAGACCCCGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCTCTCCTGGACAGACCTCGTGCAGGCGCCTCTCCTGGACCGACCTCGTGCAGGCGCCGCGCTGGACCGACCTCGTGCAGGCGCCGCGCTGGGCCATGGGGAGAGCGAGAGCCTGGTGTGCCCCTCAGGGAC".getBytes(),
                Arrays.asList(intervalOne, intervalTwo), true);
        data.add(new Object[]{contig, Arrays.asList(intervalTwo), Arrays.asList(intervalOne), 3, 1});

        intervalOne = new AlignmentInterval(new SimpleInterval("chr21", 30374719, 30375721),
                1, 1002,  TextCigarCodec.decode("966M1D36M2362H"), true, 60, 6, 960, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chr21", 30375922, 30378473),
                826, 3364,  TextCigarCodec.decode("825S33M1D1047M7D553M5D906M"), true, 60, 24, 2423, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval intervalThree = new AlignmentInterval(new SimpleInterval("chr1_KI270760v1_alt", 22529, 23531),
                1, 1002,  TextCigarCodec.decode("966M1D36M2362H"), true, 14, 3, 975, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval intervalFour = new AlignmentInterval(new SimpleInterval("chr1_KI270760v1_alt", 23681, 26220),
                826, 3364,  TextCigarCodec.decode("825H33M1D2506M"), true, 60, 2, 2517, ContigAlignmentsModifier.AlnModType.NONE);
        contig = new AlignedContig("asm027070:tig00000", "GAGCCCATCTCCTTGACTGTGGCTCTGATGCTGCCTCCACACTGGGATCTCTCTGCTCTCTTCACCTCATACCTCCTTCCCCCCACCTCACCCCATCGCCCCCGTTCTTGATCCTGCAATTGTAGAAACAGAAAGTTGGCTGATTTCTTGGGCCCGCAAATTGCCCAACAGGGAGACTGGGTGGGCGGCCCCCGCTTCCACTCCATCGCCCACCCTGATGCATCGTCTGACACTTTCAATTTATTTTTCAATTCCTCTACCATCAGAAATGACGATTAGATTTCCAGCATAAATACCGCCTTACCAAACTGAATTAATCACGGCAAGGAGGGGCACACACAGGCTCCAGCAGCCTGGGCAGAACATCCCCAGCATTAACCCTTCCGTCCTCACCCAGGCCCCCACCAGCAGGACGGAGGCTCCAGGCCTCACAGAAGACGCCACTCAAAATATCACTGGGGTCACCTAATCCCATCCCCCTTACCCTTTGCAGCCTCCCTCCTGTGGGAGTTCCTAGGAAGTGTCTTGCCCAAAGCCATCCACTCCATCAGGGCAGAGTCAGAGACACTGGCCCCTCATCTCCAGCCCCATCAGGGAAGGAGGCTCCATCCACATCCAGGACAAGATGTGGGAGTATCCGGGGTTTGGCGTTGTCCAGGACACATACGGGACGGGACTCCTGCAGACCCGAGGGTGGGGGCACCCAGTGATCACAGGGCCTGAACTGAAAGGGGTCTTGGAGAGACCTGGAGGCAGGTTCCAACCCTTGCCCCACAAACAAGACCATCACCCCTCTTTGCTGAGACTGTTCATTGCTCAGTCCAACAACCACAGCTCAGGTTGACCTCCAGCCTCCCCACTTCTCCACCTCCCTGACTCCAACCACAGCTCAGGGTGACCTCCAGCCTCCCCACTTCTCCACCTCCCTGACTCCAACCACAGCTCAGGGTGACATCCAGCCTCCCCACTTCTCCACCTCCCTGACTCCAGCCACAGCTCAGGCTCCTTCCTATGAGACCCCCATGGCCTCTCACAGCCTCTCCACTTCTATGCCTGTTCTCACCCAATCCCCATCCCTCAGCAGTCATCACCTCAAAATGCAAACACTGTCCTATGGTTTCCTGGCTCAGAACCCATCGGGCCCTCCTCTGCTCTCAAATCAGGCCCCCACCCTTCAAGGCCATGAGGACTGGGCTGGCCTGGCCCCTACCGGTCAGTGCACTCCCCCATCCTGGCTGGGTTGTCTCCTCTTTCTCCTTCAAGTTTTTCTATTTAAAATTCCCCTCCTCAGAGAACCTTCTCTGGCCACCATCCCCCAATCTAAATTAGGTTCTCCCTCCTAAGGTTCTTTCTCAAATCCATTTCCTTTCCTTCTGAGCACTTAAGCGAGCGATAATTACACACTAACTTGTGTAATTTGTTTAATAGGATCTTTGGGACAGAGACTTTATCTGACTCGCTTGATGCTGCAGCTGCTAGAACCCAGACCGTAATGTAGTGGGAGCTCAGTGCAGACTTTTGAAGGAGTAAGTGAGTAAAAGAACAACAAGCCCCTCTTGGTGCCCACCAAGTGCCAAGCTGAGACTGGGCCCTGGAGCTGGAGTCAAGATGTGGACCTGGCCTTGGTGTGCTGGGCCCTAACAGATGAGTAGGAGTTTGCCGAGCACTGAAGGTGGGGTTGACATGACCAACTTCTGAGAGGCACTCTTTGCCTCTGGATGGCCCCTTCCCAGTCACCCCAAAAGGAAGCCCTTGCCCTTTCAAAAGTGGTGAATGTGGTGGTTCAGATCGGTAGGTGTTCCTATGAATAGGTGAGGGGCCAGGCTTCAGGTCAGTTGAACCTGGGTTTGAATCCTGATTTTGCTCTTGGTACTAGGGCAGGTCACTGAGACGCTCTGAGCCTCTCTGCTCCAGGATGAGGATCCCTTCATCCATGCTCACTCAAAGTCCTGCCCACCAGGATGGAGGCAGACAGGCTGCAATGCCCTCCCCTCTCAGTGGGGGAAAAATACCAGGTCAGGCAGCCAGCAGCCGAGAATGCCAGGCAGAGCAAAGGTGTCCTAAGGGATGGACAGAATAAGGGCTTGAGAGCCTAGCCAAGGGTGAGGCTAGGAGAGGCTTCCCGGAGGACGAGGCAAGTCAGAGCTCTTTGCCTCTTACTCCCATGACTGTGGGTGCCTTTCTCCTCCTCCTCTCATTCTCTCTCCTTTCCAGCTCCTGCTCTGCTCATTTCTTCACCTCAGTCTCTCTGCCCCGACAGGAGCCCTGAGGGACACAACCCCGTCCCGAGGAATGTATCTGCCCACTTCCAGCAGGTTCCTGGAGGCCCTCTAAATTCCCCTTCCCCCCAAAGTCATCTCCCAACACTGCTGCTCCCAGGGTGGGACGCCTGCTGCTGCACCTCCACACACGTGCACACACCCAGCCAGGTGCAGACAGCGTGGGCAGTGCAGAGGGGAGGGCTGGGGATTAAGGAGTTCGTGTTCTTGAGCAGCCTGGAAAGCAGCAGGGCTTCCACAGGAGCCGCCCCTGCCCTCACCCCTGCCCAGTAGGGTTAAGGGGCTGGCTTAGATGTCACCCCAAGCCAAGGCTGTCCTTCTCAGAGGCTCCTTCCCAGCTCCCCTGAGTGGGTCAGTCCCTTCCCCTCTCTGAGCCCCTCTTTCCTCTTCTGTAAAGCAGACTCAGTGATGTTGCTCAGAGGATTGAAGGACAAAGAAAAGCAACACAATGGACAGCAGGGATTTGCAAACAGCCGGGTGCTGTACCCAAGACAGGGTATTGCTGGTGATGTCTGATGGATGGGGAGTTGAAAGACTCAGCTGTCACTGGGCAGCTGGGTCTGGTTCCCCTGAGTCATTCGTAATTCACCAACCCAGTCTATAGAAGCTTATTAAGCACTTATTGTGTGCCATGCTCCATGCAAGGGCCAAAGACACCATGAGCAGAGCCAGACCCCACCCTCAGGTTCCCCCATGGGATGGGGTTAGCCAGATGACCTGAAGGCCTCTCCAGCCAGCTCAACCCCCTTAATCCAGAATTACTCCCTGTGCCAGGCTGACGGTGTGGCCAGAGAGGCCAGGGCCTGGGAGGGGGCCTGGCAGTGGGTGGTGGGAAGAGATGGAGTGGCTGTGTCAGGGGAAGGAGAGAGCAGGTTGTTCCTGTACAGGTTTCGCTCCTCGGATAGGGGGCTGCAATGACAGCTTCCAGGAAAGACCAGGCAAGTGCCTCACCCCATCCATTCTTGCTCACCCCTGCGGCCTCTTGGCCAATGGCTGCTGTGACCCTGTCCTCCTCTGGGAATCTGGTCTCGGGGAGGAGCCCTGGACCCTGACATTGACTAGAAACCTGACCCCATGTCTGAGCA".getBytes(),
                Arrays.asList(intervalOne, intervalTwo, intervalThree, intervalFour), false);
        data.add(new Object[]{contig, Arrays.asList(intervalOne, intervalFour), Arrays.asList(intervalThree, intervalFour), 2, 2});

        // this is a case where {intervalOne} is equally good with {intervalOne, intervalTwo}, but somehow the score for latter case is tiny bit better than the first
        intervalOne = new AlignmentInterval(new SimpleInterval("chr20", 60230348, 60231029),
                1, 682, TextCigarCodec.decode("682M"), false, 57, 68, 342, ContigAlignmentsModifier.AlnModType.NONE);
        intervalTwo = new AlignmentInterval(new SimpleInterval("chrUn_JTFH01001804v1_decoy", 3674, 4300),
                1, 627, TextCigarCodec.decode("627M55H"), true, 60, 1, 622, ContigAlignmentsModifier.AlnModType.NONE);
        contig = new AlignedContig("asm005003:tig00056", "AAAACTGCTCTATCAGAAGAAAGGTTAAGCTCTGAGAGTTGAACGCACACATCACAAAGTAGTTTCTAAGAATCATTCTGTCTGGTTTTCCTATGAAGATATTGCCTTTTCTACCATAGGCCTCAAACGGCACTAAATATCCTCTTTGAAATCCTTCAAAAAGAGACTCTCAAAACTTCTCTATCGAAAGGAAGGTTCAACACCGTGAGTTGAAAGCACACATCAGAAAGAAGTTTCTGAGAAGTATTCTGTCTAGTTTTATAGGAAGAAATCACGTTTCAAAAGAAGGCCACAAAGAGGTCCAAATATCCACTTGCAGATTCTACAAAAAGAGTGTTTCAAAACTGCTCTATCAAGAGAAATGTTCATCTCCGTGAGGTGAATGCAAATATTTCAATGTAGTTTCTGACAGTGCTTCTGTCTAGTTTTTATGTGAAGATATTTCCTTTTCTACCGTAGGCCTCAAAACACTCTCAATATACACTTGCAAATTCCACAAAAAGAGTGATTCAAAACTGCTCTATCAAAAGAAATTTTAAACGCTGTAAGCTGAATGCACACATCACAAAGTAGTTTCTGAGAATGATTCTGTCTAGTTTTTCTATGAAGATATTTCCTTTTCTACCATAGGCCTTGAAGCGCTCTAAATATCCACTTGGAAATTCTACAAAAAGAGTATTTC".getBytes(),
                Arrays.asList(intervalOne, intervalTwo), true);
        data.add(new Object[]{contig, Arrays.asList(intervalOne, intervalTwo), Arrays.asList(intervalOne), 2, 1});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "contigAlignmentsHeuristicFilter", groups = "sv")
    public void testSuite(final AlignedContig contig,
                          final List<AlignmentInterval> configuration,
                          final List<AlignmentInterval> configurationEquallyGoodOrBetter,
                          final int expectedConfigurationCount,
                          final int expectedAICount) {

        final double scoreOne = AssemblyContigAlignmentsConfigPicker.computeScoreOfConfiguration(configuration, b38_canonicalChromosomes, 60);
        final double equallyGoodOrBetterScore = AssemblyContigAlignmentsConfigPicker.computeScoreOfConfiguration(configurationEquallyGoodOrBetter, b38_canonicalChromosomes, 60);
        assertTrue( scoreOne < equallyGoodOrBetterScore || scoreOne - equallyGoodOrBetterScore <= Math.ulp(equallyGoodOrBetterScore));

        assertEquals(AssemblyContigAlignmentsConfigPicker.pickBestConfigurations(contig, b38_canonicalChromosomes, 0.0).size(), expectedConfigurationCount);

        if (expectedConfigurationCount == 1) {

            final List<AlignmentInterval> alignments = AssemblyContigAlignmentsConfigPicker.gatherBestConfigurationsForOneContig(
                    SparkContextFactory.getTestSparkContext().parallelize(Collections.singletonList(contig))
                    , null, b38_seqDict, 0.0
            ).values().collect().get(0).get(0).goodMappings;
            assertEquals(alignments.size(), expectedAICount,
                    alignments.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()).toString());
        }
    }


    @DataProvider(name = "gapSplitFineTuning")
    private Object[][] createTestDataForGapSplit() {
        final List<Object[]> data = new ArrayList<>(20);

        final AlignmentInterval alignmentOne = new AlignmentInterval(new SimpleInterval("chrUn_JTFH01000492v1_decoy", 501, 1597),
                1, 1097, TextCigarCodec.decode("1097M6H"), true, 60, 1, 1092, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval alignmentTwo = new AlignmentInterval(new SimpleInterval("chr17", 26962248, 26962806),
                483, 1103, CigarUtils.invertCigar(TextCigarCodec.decode("121M1D142M1I165M62I130M482S")), false, 60, 97, 281, ContigAlignmentsModifier.AlnModType.NONE);

        final Iterable<AlignmentInterval> split = ContigAlignmentsModifier.splitGappedAlignment(alignmentTwo, GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY, 1103);
        data.add(new Object[]{
                new AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings(Arrays.asList(alignmentOne, alignmentTwo), Collections.emptyList()),
                new AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings(Collections.singletonList(alignmentOne), Lists.newArrayList(split))
        });

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "gapSplitFineTuning", groups = "sv")
    public void testGapSplit(final AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings inputConfiguration,
                             final AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings expectedOutputConfiguration) {

        final AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings configuration = AssemblyContigAlignmentsConfigPicker.splitGaps(inputConfiguration);
        Assert.assertEquals(configuration, expectedOutputConfiguration);
    }

    @Test(groups = "sv")
    public void testBreakTieByPreferringLessAlignments() {

        AlignmentInterval intervalOne = new AlignmentInterval(
                new SimpleInterval("chr21", 100000, 100100),
                1, 100, TextCigarCodec.decode("100M220S"),
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval intervalTwo = new AlignmentInterval(
                new SimpleInterval("chr21", 100099, 100122),
                99, 122, TextCigarCodec.decode("98S24M78S"),
                true, 10, 3, 241, ContigAlignmentsModifier.AlnModType.NONE);
        AlignmentInterval intervalThree = new AlignmentInterval(
                new SimpleInterval("chr21", 100121, 100200),
                122, 200,  TextCigarCodec.decode("222S78M"),
                true, 60, 0, 78, ContigAlignmentsModifier.AlnModType.NONE);

        final AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings rep1 =
                new AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings(Arrays.asList(intervalOne, intervalThree),
                        Collections.singletonList(intervalThree));
        final AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings rep2 =
                new AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings(Arrays.asList(intervalOne, intervalTwo, intervalThree),
                        Collections.emptyList());

        int threshold = 0;
        List<AssemblyContigAlignmentsConfigPicker.GoodAndBadMappings> result = AssemblyContigAlignmentsConfigPicker.breakTieByPreferringLessAlignments(Arrays.asList(rep1, rep2), threshold);
        Assert.assertEquals(IterableUtils.size(result), 2);

        threshold = 10;
        result = AssemblyContigAlignmentsConfigPicker.breakTieByPreferringLessAlignments(Arrays.asList(rep1, rep2), threshold);
        Assert.assertEquals(IterableUtils.size(result), 1);
    }

    @DataProvider(name = "forTestingNotDiscardForBadMQ")
    private Object[][] forTestingNotDiscardForBadMQ() {

        final List<Object[]> data = new ArrayList<>(20);

        final AlignedContig outForEmptyAlignments = new AlignedContig("unmapped", "AAAAAAAAAA".getBytes(), Collections.emptyList(), false);
        data.add(new Object[]{outForEmptyAlignments, false});

        final AlignedContig outForSingleBadMapping = new AlignedContig("badMap", makeDummySequence(151, (byte)'T'),
                Collections.singletonList(new AlignmentInterval(new SimpleInterval("chr1", 1000000, 1000150),
                        1, 151, TextCigarCodec.decode("151M"), true, 5, 39, 100, ContigAlignmentsModifier.AlnModType.NONE)),
                false);
        data.add(new Object[]{outForSingleBadMapping, false});

        final AlignmentInterval intervalOne = new AlignmentInterval(
                new SimpleInterval("chr21", 100000, 100100),
                1, 100, TextCigarCodec.decode("100M220S"),
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval intervalTwo = new AlignmentInterval(
                new SimpleInterval("chr21", 100099, 100122),
                99, 122, TextCigarCodec.decode("98S24M78S"),
                true, 10, 3, 241, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignedContig outForOnlyOneGaplessGoodMapping = new AlignedContig("lonelyGood", makeDummySequence(320, (byte)'C'),
                Arrays.asList(intervalOne, intervalTwo), false);
        data.add(new Object[]{outForOnlyOneGaplessGoodMapping, false});

        final AlignedContig stayForOneGappedGoodMapping =
                fromPrimarySAMRecordString("asm000576:tig00001\t16\tchr1\t30894493\t60\t106S126M67D873M\t*\t0\t0\tGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTGAGTATGAGTGGGTGAGAGTGTGAGTGGGTGAGAGTGTGTGTGGATGAGTGTGTGTGGGTGAGAGTGTGTGGGTGTATGTGTGGGCTACGTGTGTGTGGGTGTGTGTGGATGAGTGTGAGTGTGGGTGAGAGTGTGTGGGTGTGTGTGGGGAGTGTGTGTGTGTATGTGTGAATGTGTGGGCGTACACAGGCATGCACATGCACACATCCAAGCACATTCAGGGGTATTTCTGGGCGTGTGTGTCTGTCTAGGTGCATCCATGTGTGTATGCCTAGGGCTGTCTCGAGTGGCTTAGATGCCTGTGTGGGTCTGAGTGTGTCCATGTTAGCGTGGAGGGCTGTGTCTCTGAGTGGTACCTGTGCACATGGCAGTGTATGTGTGTCCCTGTCTCTGTCTCTACACGTGTGCACAGGACAAGGGTATTTGCACAGGTCCATGTTCCTCCATTTGTAAGCATGCATCTGAGCCTCCCAAGGTAGAGAGCATCTCCCAACTTCCTTAATCCACAGTGTTGGGGGTAAGCCCCTGGAGCTGAGGCAGGGGGAGGGACAGCAGATGGGGGACTCCTGATCCCCCTGGCATAGGCTACAAACACTTTTTCCTACAAAGACACCAGGGCCCCATGCATCAGCACAGACACACAACAACATGCCATCACACATGTTGCAGGTACACGCAAACCCCTAACCCAAACTAACGTCCTAAACAGAAACATGGCCTTAGCTCCACGTGGAAAAAGGATTACCTGGCCTCACTTGCAGCCTTGCAGACCATCCCTCAAGACACAGCTGAGCTCCCTCTGGATCTGAGATTGTCCCCACAAGCCCCCCAGCTCCGCAGCTCACATGAGAATAGGATCCAGCCAGCCAGGGCTCAGCTCCCAGAGGCCATGTAGGGGCTTCAGCCCCGGGGAAGCTGTTCCACTTCTCACTGCTGAGTGACCCTCAACAGCCACACTTCATCTCTGGGTTCAGTCTCCTCATCTGAACAG\t*\tSA:Z:chr2,32916450,-,107M998S,0,2;\tMD:Z:46C79^GTGTGGGTGAGAGTGTGTGGGTGTGTGTGGGGTGTGTGTGGGTGAGTGTGTGTGGATGACTGTGAGT396C40T435\tRG:Z:GATKSVContigAlignments\tNM:i:70\tAS:i:901\tXS:i:0",
                        true);
        data.add(new Object[]{stayForOneGappedGoodMapping, true});

        final AlignedContig stayForOneGoodMapping =
                fromPrimarySAMRecordString("asm022620:tig00006\t0\tchr1\t176967755\t60\t1022M\t*\t0\t0\tCTGTTGCTAAGTAGTTTAGGTCAAGGAGACCTCCTGCTCTTGAGGGGCCTGAAAATTAAGAAAATCAACTCAAATTGGAAGGATGTAGAGGAGTGCAGCGTAGACCCCACCCTGATCACATTTAGAGGCAGGAGTGCAGTGTAGACCCCACCCTGATCACATTTAGAGGCAGGAGTGTAGTGTAGACCCCATCCTGATCACATTTAGAGGCAGGAGTGTAGTGTAGACCCCATCCTGATCACATTTAAAGGCAGCAGGGCACAGTGGTTTCAGGACATGGATTCTGGCGTTAGCCTGCCTATTTCAAATGCCAGCTTAGCCACCAATCACTGTGTGACTTTGGGCAAGTTCCTTCATCATTCTGTAACCCCATATCCTCATCTGTAAAATGAGGCAGATGGTAATAACTGTATTTATTTGACAAAGTTGTTATGAGGAGGAACTAGCTAATGTAATTAAAGTGCCTAGAACAGTGCGTGGCATTTAGTGAGTGCCACACATAGAATCATTTGTTAAATAAGTTAAATTGATCAATTACATGCATCAAATATAAATACTATAAGGAGTGTTTCCCAAATTTCAATAGAGGAAGATGTGTCTTGCCCTAAGGTGCAAATGCAATGCCAACGCCAGGCTTTAACTCAGATTATAGATGTCTTCACATAGGCTGGAGAGGAGGCAACATGAGTATGGAAAGGCCAAAGCCCCGACACCCCACCTCCCACCTCCCACTTCCCACCCAGTGGTCACTTCACTGGGCAACTCACTCTACAGCATTGACATTAAACACAAAGAGCAGGGAAGTGCATATTGTGAAATGGCAAGAGAAGTCAGCGAGTAGCAAAATTGCACCAGCAAAGGAATAAAGGGGAATTTACCTAATTTCTTCCCTCTAGTCCGGAATGATCTGTCAGAGAGATGAGAAATCTTACTTGATTTATTTAGCAAATATTTGTGGAGGCCTATTGTGTGTGAACTGTGATATTAGGCAACAGGAATTGGGATAAAGAGATTCCTGCCCT\t*\tMD:Z:287A734\tRG:Z:GATKSVContigAlignments\tNM:i:1\tAS:i:1017\tXS:i:165",
                        false);
        data.add(new Object[]{stayForOneGoodMapping, true});

        final AlignedContig stayForTwoGoodMappings =
                fromPrimarySAMRecordString("asm017085:tig00001\t0\tchr10\t91874806\t60\t142S150M\t*\t0\t0\tAAAAAAAAGAATCTCTGAGAAAATACACAGGAAACTGATAATCTGGTTGCTTCTATGAAGGAAACGGCATGGCAGGGGATGGGATGGAAGGGAGGCTTGCTTTTCACTCTATTGCCTTTTGAATTCTGTACTGTGTGCATTTCACAATGTTGTGCTTAGAATTATAGAAAAACAAAAGCACAGCTCTGAACATGGGTGTCCCAATGGCACCTCAGCCTTAACCAACTTCTGCCTCTATTGCAAGCTCTTAATTGGCCTCACCTCTCACTTTCCTCATACAATTGACATTCCA\t*\tSA:Z:chr10,91873437,+,145M147S,60,0;\tMD:Z:150\tRG:Z:GATKSVContigAlignments\tNM:i:0\tAS:i:150\tXS:i:19",
                        true);
        data.add(new Object[]{stayForTwoGoodMappings, true});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forTestingNotDiscardForBadMQ", groups = "sv")
    public void testNotDiscardForBadMQ(final AlignedContig contig, final boolean shouldKeep) {
        Assert.assertEquals(AssemblyContigAlignmentsConfigPicker.notDiscardForBadMQ(contig), shouldKeep);
    }
}
