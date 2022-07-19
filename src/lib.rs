#![warn(missing_docs)]

//! Read and write minimap2 PAF files.
//! Lines are internally stored as struct [PAFLine].

/// The data structures storing PAF lines.
pub mod data;
/// Read PAF lines.
pub mod input;
/// Output PAF lines.
pub mod output;
