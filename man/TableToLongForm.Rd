\name{TableToLongForm}
\alias{TableToLongForm}
\title{
  Convert a Table to a LongForm data.frame
}
\description{
  Run TableToLongForm on the specified matrix \code{Table}.\cr
  All other arguments are optional and primarily for diagnostic
  purposes.
}
\usage{
TableToLongForm(Table, IdentResult = NULL,
                fulloutput = FALSE, diagnostics = FALSE)
}
\arguments{
  \item{Table}{
    the Table to convert, given as a character matrix.
  }
  \item{IdentResult}{
    an optional list specifying the locations of the various elements of
    the Table. By default this is automatically generated but it can be
    specified manually where the automatic detection fails.
  }
  \item{fulloutput}{
    if TRUE, returns a list containing additional information primarily
    useful for diagnostic purposes. Otherwise, and by default, the
    function only returns the converted data.frame object.
  }
  \item{diagnostics}{
    a character vector specifying the name of the file diagnostic output
    will be written to. Can also be TRUE, in which case the file name
    will be the name of the object specified in \code{Table}.
  }
}
\details{
  For details on what TableToLongForm does and what sorts of Tables it
  can convert, refer to the Technical Report found at
  \url{https://www.stat.auckland.ac.nz/~joh024/Research/TableToLongForm/TableToLongForm_TechReport.pdf}.
}
\value{
  The converted Table as a data.frame object.
}
\examples{
  ## load Toy Examples data
  data(TCData)

  ## Convert ToyExComplete
  TableToLongForm(TCData$ToyExComplete)
}