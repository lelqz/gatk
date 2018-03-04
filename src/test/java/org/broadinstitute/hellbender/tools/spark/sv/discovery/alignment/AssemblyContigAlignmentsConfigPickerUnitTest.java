package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

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
            final AlignedContig tig =
                    AssemblyContigAlignmentsConfigPicker.filterAndSplitGappedAlignmentInterval(
                            SparkContextFactory.getTestSparkContext().parallelize(Collections.singletonList(contig)),
                            null, b38_seqDict, 0.0).collect().get(0);
            assertEquals(tig.alignmentIntervals.size(), expectedAICount,
                    tig.alignmentIntervals.stream().map(AlignmentInterval::toPackedString).collect(Collectors.toList()).toString());
        }

    }


    @DataProvider(name = "gapSplitFineTuning")
    private Object[][] createTestDataForGapSplit() {
        final List<Object[]> data = new ArrayList<>(20);

        final AlignmentInterval alignmentOne = new AlignmentInterval(new SimpleInterval("chrUn_JTFH01000492v1_decoy", 501, 1597),
                1, 1097, TextCigarCodec.decode("1097M6H"), true, 60, 1, 1092, ContigAlignmentsModifier.AlnModType.NONE);
        final AlignmentInterval alignmentTwo = new AlignmentInterval(new SimpleInterval("chr17", 26962248, 26962806),
                483, 1103, CigarUtils.invertCigar(TextCigarCodec.decode("121M1D142M1I165M62I130M482S")), false, 60, 97, 281, ContigAlignmentsModifier.AlnModType.NONE);

        data.add(new Object[]{Arrays.asList(alignmentOne, alignmentTwo), Collections.singletonList(alignmentOne)});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "gapSplitFineTuning", groups = "sv")
    public void testGapSplit(final List<AlignmentInterval> inputConfiguration, final List<AlignmentInterval> expectedOutputConfiguration) {

        final List<AlignmentInterval> alignmentIntervals = AssemblyContigAlignmentsConfigPicker.splitGaps(inputConfiguration);
        Assert.assertEquals(alignmentIntervals, expectedOutputConfiguration);
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
}
