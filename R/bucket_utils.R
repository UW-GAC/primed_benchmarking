#' Copy a file from a Google Cloud Storage bucket to a local path
#'
#' Downloads a file from a Google Cloud Storage (GCS) bucket to a local path
#' using \code{gsutil cp}. If the destination file already exists locally, the
#' copy is skipped and a message is printed instead.
#'
#' @param gspath Character. The GCS path of the file to copy
#'   (e.g. \code{"gs://my-bucket/path/to/file.txt"}).
#' @param newfilename Character. The local destination path where the file
#'   should be saved.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @details
#' This function requires that \code{gsutil} (part of the Google Cloud SDK)
#' is installed and available on the system \code{PATH}. An error is thrown if
#' \code{gsutil} is not found or if the \code{gsutil cp} command exits with a
#' non-zero status. Stdout and stderr from the command are captured and printed
#' to assist with debugging.
#'
#' @examples
#' \dontrun{
#' get_file_from_bucket(
#'   gspath      = "gs://my-bucket/data/scores.txt",
#'   newfilename = "scores.txt"
#' )
#' }
#'
#' @export
get_file_from_bucket <- function(gspath, newfilename) {
    if (!is.character(gspath) || length(gspath) != 1 || nchar(gspath) == 0) {
        stop("'gspath' must be a non-empty character string.")
    }
    if (!is.character(newfilename) || length(newfilename) != 1 ||
            nchar(newfilename) == 0) {
        stop("'newfilename' must be a non-empty character string.")
    }

    .check_gsutil()

    if (!file.exists(newfilename)) {
        output <- system(
            paste0("gsutil cp ", gspath, " ", newfilename, " 2>&1"),
            intern = TRUE
        )
        print(output)
        status <- attr(output, "status")
        if (!is.null(status) && status != 0) {
            stop(sprintf(
                "gsutil cp failed with exit status %d. ",
                status
            ))
        }
    } else {
        print(paste(newfilename, "already copied"))
    }
    invisible(NULL)
}


#' Copy a local file to a Google Cloud Storage bucket
#'
#' Uploads a local file to a Google Cloud Storage (GCS) bucket using
#' \code{gsutil cp}, then lists the uploaded file with \code{gsutil ls -l}
#' to confirm the transfer.
#'
#' @param filename Character. Path to the local file to upload.
#' @param gspath Character. The GCS bucket or prefix path to copy the file
#'   into (e.g. \code{"gs://my-bucket/path"}).
#' @param newfilename Character. The name to give the file in the bucket.
#'   The file will be written to \code{<gspath>/<newfilename>}.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @details
#' This function requires that \code{gsutil} (part of the Google Cloud SDK)
#' is installed and available on the system \code{PATH}. An error is thrown if
#' \code{gsutil} is not found, if the local \code{filename} does not exist, or
#' if either \code{gsutil} command exits with a non-zero status. Stdout and
#' stderr from each command are captured and printed to assist with debugging.
#'
#' @examples
#' \dontrun{
#' copy_file_to_bucket(
#'   filename    = "results.txt",
#'   gspath      = "gs://my-bucket/output",
#'   newfilename = "results.txt"
#' )
#' }
#'
#' @export
copy_file_to_bucket <- function(filename, gspath, newfilename) {
    if (!is.character(filename) || length(filename) != 1 ||
            nchar(filename) == 0) {
        stop("'filename' must be a non-empty character string.")
    }
    if (!is.character(gspath) || length(gspath) != 1 || nchar(gspath) == 0) {
        stop("'gspath' must be a non-empty character string.")
    }
    if (!is.character(newfilename) || length(newfilename) != 1 ||
            nchar(newfilename) == 0) {
        stop("'newfilename' must be a non-empty character string.")
    }

    if (!file.exists(filename)) {
        stop(sprintf("File '%s' does not exist.", filename))
    }

    .check_gsutil()

    dest <- paste0(gspath, "/", newfilename)

    output_cp <- system(
        paste0("gsutil cp ", filename, " ", dest, " 2>&1"),
        intern = TRUE
    )
    print(output_cp)
    status_cp <- attr(output_cp, "status")
    if (!is.null(status_cp) && status_cp != 0) {
        stop(sprintf(
            "gsutil cp failed with exit status %d. ",
            status_cp
        ))
    }

    output_ls <- system(
        paste0("gsutil ls -l ", dest, " 2>&1"),
        intern = TRUE
    )
    print(output_ls)
    status_ls <- attr(output_ls, "status")
    if (!is.null(status_ls) && status_ls != 0) {
        stop(sprintf(
            "gsutil ls failed with exit status %d. ",
            status_ls
        ))
    }

    invisible(NULL)
}


# ---- Internal helpers --------------------------------------------------------

#' @keywords internal
.check_gsutil <- function() {
    path <- Sys.which("gsutil")
    if (nchar(path) == 0) {
        stop(
            "gsutil is not installed or not found on PATH. ",
            "Please install the Google Cloud SDK: ",
            "https://cloud.google.com/sdk/docs/install"
        )
    }
}
