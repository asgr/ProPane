\name{propaneFrameFinder}
\alias{propaneFrameFinder}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Finder Frame Coverage
}
\description{
For a given region of sky and a number of on-disk FITS files, this function will find all nearby frames and optional plot their coverage. This is useful when trying to figure out an optimal WCS for stacking purposes.
}
\usage{
propaneFrameFinder(filelist = NULL, dirlist = NULL, extlist = 1, RAcen = 180, Deccen = 0,
  rad = 1, cores = 1, plot = TRUE, proj = TRUE, border = "red", col = hsv(alpha = 0.2),
  ...)
}
\arguments{
  \item{filelist}{
Character vector; vector of full paths of FITS files to analyse.
}
  \item{dirlist}{
Character vector; vector of full paths of directories to analyse.
}
  \item{extlist}{
Integer vector; the extensions to use. If length 1 then it will be used for all \option{filelist}, but otherwise it should be the same length as the final \option{filelist}
}
  \item{RAcen}{
Numeric scalar; the Right Ascension of the centre of the field in degrees.
}
  \item{Deccen}{
Numeric scalar; the Declination of the centre of the field in degrees.
}
  \item{rad}{
Numeric scalar; the radius from \option{RAcen}, \option{Deccen} to search in units of \option{radunit}.
}
  \item{cores}{
Integer scalar; the number of cores to run on.
}
  \item{plot}{
Logical; should a coverage plot be generated?
}
  \item{proj}{
If \option{plot} = TRUE then should \code{magproj} be used for more slower but more accurate projection?
}
  \item{border}{
Scalar; colour of the border of frames when \option{plot} = TRUE.
}
  \item{col}{
Scalar; colour inside of frames when \option{plot} = TRUE.
}
  \item{\dots}{
Other arguments to pass to \code{\link{Rfits_key_scan}}.
}
}
\details{
At the simplest level this function uses \code{\link{Rfits_key_scan}} to get useful information about target FITS files. It then uses \code{\link{coordmatchsing}} to find frame centres that are within the circle specified by \option{RAcen}, \option{Deccen} and \option{rad} and trims the output of \code{\link{Rfits_key_scan}} to just these frames. Optionally it can then plot the coverage.
}
\value{
The output is the trimmed the output of \code{\link{Rfits_key_scan}} pluts a separation column (sep) in units of \option{radunit}.
}
\author{
Aaron Robotham
}
\note{
It is not necessary to trim down inputs before running \code{\link{propaneStackWarpInVar}}, but this function makes it easy to find an optimal WCS for the stacked frames.
}
\seealso{
\code{\link{Rfits_key_scan}}
}
\examples{
ex_frames = propaneFrameFinder(dirlist = system.file('extdata/stack/', package="ProPane"),
  extlist = 2, RAcen=36.9, Deccen=-5.2, rad=0.1)

print(ex_frames)
}

