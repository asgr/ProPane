\name{propaneSegimShare}
\alias{propaneSegimShare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate Flux Share Matrix
}
\description{
Given a higher resolution input segmentation map, this function calculates how the original segments are split across a lower resolution segmentation map (\option{segim_warp}). The output of this will be used to run \code{profoundProFound} using a new target segmention map, and the outputs from this photometry run and the \option{sharemat} provided will then be passed to \code{profoundShareFlux} to create fluxes for sources that are not present in \option{segim_warp} but exist in the original \option{segim}.
}
\usage{
propaneSegimShare(segim = NULL, keyvalues_in = NULL, keyvalues_out = NULL, pixcut = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segim}{
Integer matrix; required, the segmentation map we want to remap. If \option{segim_in} is a list as created by \code{Rfits_read_image} then the image part of the list is parsed to \option{segim_in} and the correct header part is passed to \option{header_in}.
}
  \item{keyvalues_in}{
List; input header values to be used for the WCS.
}
  \item{keyvalues_out}{
List; output header values to be used for the WCS.
}
  \item{pixcut}{
Integer scalar; the number of pixels required to identify an object on the new low resolution segmentation map.
}
}
\value{
A list containing:

\item{segID_in }{Vector, the input segment IDs present in\option{segim_in}.}
\item{segID_warp }{Vector, the output segment IDs present in\option{segim_warp}.}
\item{segim_warp}{Integer matrix; the remapped image using the target WCS.}
\item{sharemat}{Numeric matrix; the sharing matrix which will have dimenstion length(segID_warp) x length(segID_in).}
\item{shareseg}{Numeric vector; the fraction of non-sky \option{segim_in} pixels that match the output \option{segim_warp}, so if near 1 then \option{segim_warp} segment should be quite clean of subpixel contaminants.}
}
\author{
Aaron Robotham
}
\seealso{
\code{\link{propaneSegimWarp}}
}
\examples{
#None yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ segim }% use one of  RShowDoc("KEYWORDS")
