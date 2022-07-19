/// The result type of this crate using the error type of this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// All errors that can occur while parsing.
///
/// Printing is error-free, except for IO errors.
/// Therefore when printing, Rust's [std::io::Error] is used.
pub enum Error {
    /// Parsing error when parsing the column into the expected type.
    ColumnParseError(Box<dyn std::error::Error>),

    /// The line was ended, but it was not expected (e.g. further columns are expected instead).
    UnexpectedEndOfLine,

    /// The file was ended, but it was not expected (e.g. further columns are expected instead).
    UnexpectedEndOfFile,

    /// An unexpected character was found.
    UnexpectedCharacter,

    /// An optional column with an unexpected name or type was found.
    UnexpectedOptionalColumn,

    /// A cigar string could not be parsed.
    MalformedCigar,

    /// An alignment difference string could not be parsed.
    MalformedAlignmentDifference,
}

pub(crate) fn column_parse_error<SourceError: 'static + std::error::Error>(
    source_error: SourceError,
) -> Error {
    Error::ColumnParseError(Box::new(source_error))
}
