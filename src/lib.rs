#![warn(missing_docs)]
#![allow(clippy::iter_nth_zero)]

//! Read and write minimap2 PAF files.
//! Lines are internally stored as struct [PAFLine].

/// The data structures storing PAF lines.
pub mod data;
/// Error handling types.
pub mod error;
/// Read PAF lines.
pub mod input;
/// Output PAF lines.
pub mod output;
