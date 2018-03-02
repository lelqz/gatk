package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import htsjdk.samtools.SAMFileHeader;

import java.util.List;

public class AnnotatedIntervalHeader {
    private final String contigColumnName;
    private final String startColumnName;
    private final String endColumnName;
    private final List<String> annotations;
    private final SAMFileHeader samFileHeader;
    private final List<String> comments;

    /**
     * @param samFileHeader SAM file header as a multiline string.  {@code null} is allowed, if not available.
     * @param comments Comments to prepend to the xsv file.  Use an empty list, if no comments are needed.  Never {@code null}.
     * @param annotations annotation names that do not include the locatable column names.  Never {@code null}.
     * @param contigColumnName how contig should be rendered.  Never {@code null}.
     * @param startColumnName how start position should be rendered.  Never {@code null}.
     * @param endColumnName how end position should be rendered.  Never {@code null}.
     */
    public AnnotatedIntervalHeader(final String contigColumnName, final String startColumnName, final String endColumnName, final List<String> annotations, final SAMFileHeader samFileHeader, final List<String> comments) {
        this.contigColumnName = contigColumnName;
        this.startColumnName = startColumnName;
        this.endColumnName = endColumnName;
        this.annotations = annotations;
        this.samFileHeader = samFileHeader;
        this.comments = comments;
    }


    public String getContigColumnName() {
        return contigColumnName;
    }

    public String getStartColumnName() {
        return startColumnName;
    }

    public String getEndColumnName() {
        return endColumnName;
    }

    public List<String> getAnnotations() {
        return annotations;
    }

    public SAMFileHeader getSamFileHeader() {
        return samFileHeader;
    }

    public List<String> getComments() {
        return comments;
    }
}
