/// The result type of this crate using the error type of this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// All errors that can occur while parsing.
///
/// Printing is error-free, except for IO errors.
/// Therefore when printing, Rust's [std::io::Error] is used.
#[derive(Debug)]
pub enum Error {
    /// An I/O error.
    IOError(std::io::Error),

    /// Parsing error when parsing the column into the expected type.
    ColumnParseError,

    /// The line was ended, but it was not expected (e.g. further columns are expected instead).
    UnexpectedEndOfLine,

    /// The file was ended, but it was not expected (e.g. further columns are expected instead).
    UnexpectedEndOfFile,

    /// An unexpected character was found.
    UnexpectedCharacter,

    /// An optional column with an unexpected name or type was found.
    UnexpectedOptionalColumn {
        /// The header of the column, or a part of it.
        /// The header is the name and type.
        column_header: String,
    },

    /// A cigar string could not be parsed.
    MalformedCigar,

    /// An alignment difference string could not be parsed.
    MalformedAlignmentDifference,

    /// Quick and dirty: simply use strings to report errors.
    Message(String),
}

impl From<std::io::Error> for Error {
    fn from(error: std::io::Error) -> Self {
        Self::IOError(error)
    }
}
