/// A line in a minimap PAF file.
///
/// The field names are adapted from the [minimap2 man-page][1]. Check it out for more info.
///
/// [1]: https://lh3.github.io/minimap2/minimap2.html#10
#[derive(Clone, Debug, PartialEq)]
#[allow(missing_docs)]
pub struct PAFLine {
    // required fields
    pub query_sequence_name: String,
    pub query_sequence_length: usize,
    pub query_start_coordinate: usize,
    pub query_end_coordinate: usize,
    pub strand: bool,
    pub target_sequence_name: String,
    pub target_sequence_length: usize,
    pub target_start_coordinate_on_original_strand: usize,
    pub target_end_coordinate_on_original_strand: usize,
    pub number_of_matching_bases: usize,
    pub number_of_bases_and_gaps: usize,
    pub mapping_quality: u8,

    // optional fields
    pub alignment_type: Option<AlignmentType>,
    pub number_of_minimisers: Option<usize>,
    pub chaining_score: Option<usize>,
    pub best_secondary_chaining_score: Option<usize>,
    pub total_number_of_mismatches_and_gaps: Option<usize>,
    pub unknown_md: Option<String>,
    pub dp_alignment_score: Option<usize>,
    pub supplementary_alignments: Option<String>,
    pub best_segment_dp_score: Option<usize>,
    pub number_of_ambiguous_bases: Option<usize>,
    pub transcript_strand: Option<String>,
    pub cigar_string: Option<Cigar>,
    pub difference_string: Option<AlignmentDifference>,
    pub approximate_per_base_sequence_divergence: Option<f64>,
    pub gap_compressed_per_base_sequence_divergence: Option<f64>,
    pub length_of_query_regions_with_repetitive_seeds: Option<usize>,
}

/// The type of a minimap2 alignment. See the [minimap2 readme](https://github.com/lh3/minimap2#algorithm-overview) for more information.
#[allow(missing_docs)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum AlignmentType {
    Primary,
    Secondary,
    PrimaryInversion,
    SecondaryInversion,
}

/// A CIGAR string. See the [this page](https://www.drive5.com/usearch/manual/cigar.html) for more information.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Cigar(pub Vec<CigarColumn>);

/// A column of a CIGAR string. See the [this page](https://www.drive5.com/usearch/manual/cigar.html) for more information.
#[allow(missing_docs)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CigarColumn {
    Match(usize),
    Insertion(usize),
    Deletion(usize),
    Mismatch(usize),
}

/// An alignment difference string. See the [minimap2 man-page](https://lh3.github.io/minimap2/minimap2.html#10) for more information.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AlignmentDifference(pub Vec<DifferenceColumn>);

/// A column of a difference string. See the [minimap2 man-page](https://lh3.github.io/minimap2/minimap2.html#10) for more information.
#[allow(missing_docs)]
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DifferenceColumn {
    Match {
        length: usize,
    },
    Insertion {
        superfluous_query_characters: String,
    },
    Deletion {
        missing_query_characters: String,
    },
    Mismatch {
        reference: char,
        query: char,
    },
}
