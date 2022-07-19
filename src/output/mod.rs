use crate::data::{
    AlignmentDifference, AlignmentType, Cigar, CigarColumn, DifferenceColumn, PAFLine,
};
use std::fmt::{Display, Formatter};

impl Display for PAFLine {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        // required fields
        write!(f, "{}", self.query_sequence_name)?;
        write!(f, "\t{}", self.query_sequence_length)?;
        write!(f, "\t{}", self.query_start_coordinate)?;
        write!(f, "\t{}", self.query_end_coordinate)?;
        write!(f, "\t{}", if self.strand { '+' } else { '-' })?;
        write!(f, "\t{}", self.target_sequence_name)?;
        write!(f, "\t{}", self.target_sequence_length)?;
        write!(f, "\t{}", self.target_start_coordinate_on_original_strand)?;
        write!(f, "\t{}", self.target_end_coordinate_on_original_strand)?;
        write!(f, "\t{}", self.number_of_matching_bases)?;
        write!(f, "\t{}", self.number_of_bases_and_gaps)?;
        write!(f, "\t{}", self.mapping_quality)?;

        // optional fields with known order
        if let Some(x) = self.total_number_of_mismatches_and_gaps {
            write!(f, "\tNM:i:{x}")?;
        }
        if let Some(x) = self.best_segment_dp_score {
            write!(f, "\tms:i:{x}")?;
        }
        if let Some(x) = self.dp_alignment_score {
            write!(f, "\tAS:i:{x}")?;
        }
        if let Some(x) = self.number_of_ambiguous_bases {
            write!(f, "\tnn:i:{x}")?;
        }
        if let Some(x) = &self.alignment_type {
            write!(f, "\ttp:A:{x}")?;
        }
        if let Some(x) = self.number_of_minimisers {
            write!(f, "\tcm:i:{x}")?;
        }
        if let Some(x) = self.chaining_score {
            write!(f, "\ts1:i:{x}")?;
        }
        if let Some(x) = self.best_secondary_chaining_score {
            write!(f, "\ts2:i:{x}")?;
        }
        if let Some(x) = self.gap_compressed_per_base_sequence_divergence {
            write!(f, "\tde:f:{x}")?;
        }
        if let Some(x) = self.length_of_query_regions_with_repetitive_seeds {
            write!(f, "\trl:i:{x}")?;
        }

        // optional fields without known order, or that do not contradict the known order (cg and cs)
        // the order is the same as in the man page at https://lh3.github.io/minimap2/minimap2.html#10
        if let Some(x) = &self.unknown_md {
            write!(f, "\tMD:Z:{x}")?;
        }
        if let Some(x) = &self.supplementary_alignments {
            write!(f, "\tSA:Z:{x}")?;
        }
        if let Some(x) = &self.transcript_strand {
            write!(f, "\tts:A:{x}")?;
        }
        if let Some(x) = &self.cigar_string {
            write!(f, "\tcg:Z:{x}")?;
        }
        if let Some(x) = &self.difference_string {
            write!(f, "\tcs:Z:{x}")?;
        }
        if let Some(x) = self.approximate_per_base_sequence_divergence {
            write!(f, "\tdv:f:{x}")?;
        }

        Ok(())
    }
}

impl Display for AlignmentType {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                AlignmentType::Primary => "P",
                AlignmentType::Secondary => "S",
                AlignmentType::PrimaryInversion => "I",
                AlignmentType::SecondaryInversion => "i",
            }
        )
    }
}

impl Display for Cigar {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for cigar_column in &self.0 {
            write!(
                f,
                "{}{}",
                match cigar_column {
                    CigarColumn::Match(length)
                    | CigarColumn::Insertion(length)
                    | CigarColumn::Deletion(length)
                    | CigarColumn::Mismatch(length) => length,
                },
                match cigar_column {
                    CigarColumn::Match(_) => "M",
                    CigarColumn::Insertion(_) => "I",
                    CigarColumn::Deletion(_) => "D",
                    CigarColumn::Mismatch(_) => "X",
                }
            )?;
        }

        Ok(())
    }
}

impl Display for AlignmentDifference {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for difference_column in &self.0 {
            match difference_column {
                DifferenceColumn::Match { length } => write!(f, ":{length}")?,
                DifferenceColumn::Insertion {
                    superfluous_query_characters,
                } => write!(f, "+{superfluous_query_characters}")?,
                DifferenceColumn::Deletion {
                    missing_query_characters,
                } => write!(f, "-{missing_query_characters}")?,
                DifferenceColumn::Mismatch { reference, query } => {
                    write!(f, "*{reference}{query}")?
                }
            }
        }

        Ok(())
    }
}
