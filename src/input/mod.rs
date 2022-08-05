use crate::data::{
    AlignmentDifference, AlignmentType, Cigar, CigarColumn, DifferenceColumn, PAFLine,
};
use crate::error::{Error, Result};
use log::warn;
use std::io::{BufRead, BufReader, Read};
use std::str::FromStr;

/// Parse all lines in the given input string.
pub fn parse_input_str(mut string: &str) -> Result<Vec<PAFLine>> {
    let mut result = Vec::new();
    while !string.is_empty() {
        result.push(parse_line(&mut string)?);
    }
    Ok(result)
}

/// Parse all lines in the given input.
/// Terminate when the input has no further lines.
pub fn parse_input_read<Input: Read>(input: &mut Input) -> Result<Vec<PAFLine>> {
    let mut result = Vec::new();
    let reader = BufReader::new(input);

    for line in reader.lines() {
        result.push(parse_line(&mut line?.as_str())?);
    }

    Ok(result)
}

/// Parses a line of input into a [PAFLine].
/// The given string slice is advanced past the parsed line of input.
pub fn parse_line(string: &mut &str) -> Result<PAFLine> {
    // required fields
    let query_sequence_name = parse_column(string, false)?;
    let query_sequence_length = parse_column(string, false)?;
    let query_start_coordinate = parse_column(string, false)?;
    let query_end_coordinate = parse_column(string, false)?;
    let strand: String = parse_column(string, false)?;
    let strand = if strand == "+" {
        true
    } else if strand == "-" {
        false
    } else {
        return Err(Error::UnexpectedCharacter);
    };
    let target_sequence_name = parse_column(string, false)?;
    let target_sequence_length = parse_column(string, false)?;
    let target_start_coordinate_on_original_strand = parse_column(string, false)?;
    let target_end_coordinate_on_original_strand = parse_column(string, false)?;
    let number_of_matching_bases = parse_column(string, false)?;
    let number_of_bases_and_gaps = parse_column(string, false)?;
    let mapping_quality = parse_column(string, true)?;

    // optional fields
    let mut alignment_type = None;
    let mut number_of_minimisers = None;
    let mut chaining_score = None;
    let mut best_secondary_chaining_score = None;
    let mut total_number_of_mismatches_and_gaps = None;
    let mut unknown_md = None;
    let mut dp_alignment_score = None;
    let mut supplementary_alignments = None;
    let mut best_segment_dp_score = None;
    let mut number_of_ambiguous_bases = None;
    let mut transcript_strand = None;
    let mut cigar_string = None;
    let mut difference_string = None;
    let mut approximate_per_base_sequence_divergence = None;
    let mut gap_compressed_per_base_sequence_divergence = None;
    let mut length_of_query_regions_with_repetitive_seeds = None;
    let mut unknown_fields = Vec::new();

    loop {
        if string.starts_with('\n') {
            *string = &string[1..];
            break;
        }

        if string.len() < 6 {
            if !string.is_empty() && *string != "\n" {
                return Err(Error::UnexpectedCharacter);
            } else {
                break;
            }
        }

        let header = &string[..5];
        *string = &string[5..];
        match header {
            "tp:A:" => {
                alignment_type = Some(match extract_column_value(string)? {
                    "P" => AlignmentType::Primary,
                    "S" => AlignmentType::Secondary,
                    "I" => AlignmentType::PrimaryInversion,
                    "i" => AlignmentType::SecondaryInversion,
                    _ => return Err(Error::UnexpectedCharacter),
                });
            }
            "cm:i:" => {
                number_of_minimisers = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed cm:i:".to_string()))?,
                )
            }
            "s1:i:" => {
                chaining_score = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed s1:i:".to_string()))?,
                )
            }
            "s2:i:" => {
                best_secondary_chaining_score = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed s2:i:".to_string()))?,
                )
            }
            "NM:i:" => {
                total_number_of_mismatches_and_gaps = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed NM:i:".to_string()))?,
                )
            }
            "MD:Z:" => {
                unknown_md = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed MD:Z:".to_string()))?,
                )
            }
            "AS:i:" => {
                dp_alignment_score = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed AS:i:".to_string()))?,
                )
            }
            "SA:Z:" => {
                supplementary_alignments = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed SA:Z:".to_string()))?,
                )
            }
            "ms:i:" => {
                best_segment_dp_score = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed ms:i:".to_string()))?,
                )
            }
            "nn:i:" => {
                number_of_ambiguous_bases = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed nn:i:".to_string()))?,
                )
            }
            "ts:A:" => {
                transcript_strand = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed ts:i:".to_string()))?,
                )
            }
            "cg:Z:" => cigar_string = Some(parse_cigar(string)?),
            "cs:Z:" => difference_string = Some(parse_alignment_difference(string)?),
            "dv:f:" => {
                approximate_per_base_sequence_divergence = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed dv:f:".to_string()))?,
                )
            }
            "de:f:" => {
                gap_compressed_per_base_sequence_divergence = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed de:f:".to_string()))?,
                )
            }
            "rl:i:" => {
                length_of_query_regions_with_repetitive_seeds = Some(
                    extract_column_value(string)?
                        .parse()
                        .map_err(|_| Error::Message("Malformed rl:i:".to_string()))?,
                )
            }
            other => {
                warn!("Found unknown field: {other}");
                let mut unknown_field = other.to_string();
                unknown_field.push_str(extract_column_value(string)?);
                unknown_fields.push(unknown_field);
            }
        }
    }

    Ok(PAFLine {
        // required fields
        query_sequence_name,
        query_sequence_length,
        query_start_coordinate,
        query_end_coordinate,
        strand,
        target_sequence_name,
        target_sequence_length,
        target_start_coordinate_on_original_strand,
        target_end_coordinate_on_original_strand,
        number_of_matching_bases,
        number_of_bases_and_gaps,
        mapping_quality,

        // optional fields
        alignment_type,
        number_of_minimisers,
        chaining_score,
        best_secondary_chaining_score,
        total_number_of_mismatches_and_gaps,
        unknown_md,
        dp_alignment_score,
        supplementary_alignments,
        best_segment_dp_score,
        number_of_ambiguous_bases,
        transcript_strand,
        cigar_string,
        difference_string,
        approximate_per_base_sequence_divergence,
        gap_compressed_per_base_sequence_divergence,
        length_of_query_regions_with_repetitive_seeds,

        unknown_fields,
    })
}

fn parse_column<Type: FromStr>(string: &mut &str, allow_eol: bool) -> Result<Type>
where
    Type::Err: 'static + std::error::Error,
{
    let limit = if let Some(limit) = string.find(['\t', '\n']) {
        if string.chars().nth(limit).unwrap() == '\n' && !allow_eol {
            return Err(Error::UnexpectedEndOfLine);
        }
        limit
    } else {
        if !allow_eol {
            return Err(Error::UnexpectedEndOfFile);
        }
        string.len()
    };

    let column = &string[..limit];
    *string = &string[limit + 1..];

    column.parse().map_err(|_| Error::ColumnParseError)
}

fn extract_column_value<'input: 'output, 'output>(
    string: &mut &'input str,
) -> Result<&'output str> {
    let limit = if let Some(limit) = string.find(['\t', '\n']) {
        limit
    } else {
        string.len()
    };

    let column = &string[..limit];
    *string = &string[(limit + 1).min(string.len())..];

    Ok(column)
}

fn parse_cigar(string: &mut &str) -> Result<Cigar> {
    let mut result = Vec::new();

    while !string.is_empty() && !string.starts_with(['\t', '\n']) {
        let limit = if let Some(limit) = string.find(['M', 'D', 'I', 'X']) {
            limit
        } else {
            return Err(Error::MalformedCigar);
        };

        let count = string[..limit]
            .parse()
            .map_err(|_| Error::Message(format!("Malformed cigar count: {}", &string[..limit])))?;
        result.push(match &string[limit..limit + 1] {
            "M" => CigarColumn::Match(count),
            "D" => CigarColumn::Deletion(count),
            "I" => CigarColumn::Insertion(count),
            "X" => CigarColumn::Mismatch(count),
            _ => return Err(Error::MalformedCigar),
        });
        *string = &string[limit + 1..];
    }

    if string.starts_with(['\t']) {
        *string = &string[1..];
    }

    Ok(Cigar(result))
}

fn parse_alignment_difference(string: &mut &str) -> Result<AlignmentDifference> {
    let mut result = Vec::new();

    while !string.is_empty() && !string.starts_with(['\t', '\n']) {
        let limit = string[1..]
            .find([':', '-', '+', '*', '\t', '\n'])
            .map(|limit| limit + 1)
            .unwrap_or(string.len());
        let marker = &string[..1];
        let characters = &string[1..limit];
        *string = &string[limit..];

        result.push(match marker {
            ":" => DifferenceColumn::Match {
                length: characters.parse().map_err(|_| {
                    Error::Message(format!(
                        "Malformed alignment difference length: {}",
                        characters
                    ))
                })?,
            },
            "-" => DifferenceColumn::Deletion {
                missing_query_characters: characters.to_string(),
            },
            "+" => DifferenceColumn::Insertion {
                superfluous_query_characters: characters.to_string(),
            },
            "*" => {
                if characters.chars().count() != 2 {
                    return Err(Error::MalformedAlignmentDifference);
                }
                DifferenceColumn::Mismatch {
                    reference: characters.chars().nth(0).unwrap(),
                    query: characters.chars().nth(1).unwrap(),
                }
            }
            _ => return Err(Error::MalformedAlignmentDifference),
        })
    }

    if string.starts_with(['\t']) {
        *string = &string[1..];
    }

    Ok(AlignmentDifference(result))
}
