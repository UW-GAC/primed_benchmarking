context("Bucket utility functions")

# Tests for get_file_from_bucket -----------------------------------------------

test_that("get_file_from_bucket errors on non-character gspath", {
    expect_error(
        get_file_from_bucket(123, "file.txt"),
        "'gspath' must be a non-empty character string"
    )
})

test_that("get_file_from_bucket errors on empty gspath", {
    expect_error(
        get_file_from_bucket("", "file.txt"),
        "'gspath' must be a non-empty character string"
    )
})

test_that("get_file_from_bucket errors on length-0 gspath", {
    expect_error(
        get_file_from_bucket(character(0), "file.txt"),
        "'gspath' must be a non-empty character string"
    )
})

test_that("get_file_from_bucket errors on non-character newfilename", {
    expect_error(
        get_file_from_bucket("gs://bucket/file.txt", 123),
        "'newfilename' must be a non-empty character string"
    )
})

test_that("get_file_from_bucket errors on empty newfilename", {
    expect_error(
        get_file_from_bucket("gs://bucket/file.txt", ""),
        "'newfilename' must be a non-empty character string"
    )
})

test_that("get_file_from_bucket errors when gsutil is not installed", {
    mockery::stub(get_file_from_bucket, "Sys.which", function(x) "")
    expect_error(
        get_file_from_bucket("gs://bucket/file.txt", "file.txt"),
        "gsutil is not installed"
    )
})

test_that("get_file_from_bucket prints 'already copied' when file exists", {
    mockery::stub(get_file_from_bucket, "Sys.which",
                  function(x) "/usr/bin/gsutil")
    mockery::stub(get_file_from_bucket, "file.exists", function(x) TRUE)
    expect_output(
        get_file_from_bucket("gs://bucket/file.txt", "file.txt"),
        "already copied"
    )
})

test_that("get_file_from_bucket calls gsutil when file does not exist", {
    mockery::stub(get_file_from_bucket, "Sys.which",
                  function(x) "/usr/bin/gsutil")
    mockery::stub(get_file_from_bucket, "file.exists", function(x) FALSE)
    mockery::stub(get_file_from_bucket, "system", function(cmd, intern) {
        expect_true(grepl("gsutil cp", cmd))
        character(0)
    })
    expect_invisible(
        get_file_from_bucket("gs://bucket/file.txt", "file.txt")
    )
})

test_that("get_file_from_bucket errors when gsutil cp exits non-zero", {
    mockery::stub(get_file_from_bucket, "Sys.which",
                  function(x) "/usr/bin/gsutil")
    mockery::stub(get_file_from_bucket, "file.exists", function(x) FALSE)
    mockery::stub(get_file_from_bucket, "system", function(cmd, intern) {
        out <- "CommandException: No URLs matched"
        attr(out, "status") <- 1L
        out
    })
    expect_error(
        get_file_from_bucket("gs://bucket/file.txt", "file.txt"),
        "gsutil cp failed"
    )
})

# Tests for copy_file_to_bucket ------------------------------------------------

test_that("copy_file_to_bucket errors on non-character filename", {
    expect_error(
        copy_file_to_bucket(123, "gs://bucket", "out.txt"),
        "'filename' must be a non-empty character string"
    )
})

test_that("copy_file_to_bucket errors on empty filename", {
    expect_error(
        copy_file_to_bucket("", "gs://bucket", "out.txt"),
        "'filename' must be a non-empty character string"
    )
})

test_that("copy_file_to_bucket errors on non-character gspath", {
    expect_error(
        copy_file_to_bucket("local.txt", 123, "out.txt"),
        "'gspath' must be a non-empty character string"
    )
})

test_that("copy_file_to_bucket errors on empty gspath", {
    expect_error(
        copy_file_to_bucket("local.txt", "", "out.txt"),
        "'gspath' must be a non-empty character string"
    )
})

test_that("copy_file_to_bucket errors on non-character newfilename", {
    expect_error(
        copy_file_to_bucket("local.txt", "gs://bucket", 123),
        "'newfilename' must be a non-empty character string"
    )
})

test_that("copy_file_to_bucket errors on empty newfilename", {
    expect_error(
        copy_file_to_bucket("local.txt", "gs://bucket", ""),
        "'newfilename' must be a non-empty character string"
    )
})

test_that("copy_file_to_bucket errors when local file does not exist", {
    mockery::stub(copy_file_to_bucket, "file.exists", function(x) FALSE)
    expect_error(
        copy_file_to_bucket("missing.txt", "gs://bucket", "out.txt"),
        "does not exist"
    )
})

test_that("copy_file_to_bucket errors when gsutil is not installed", {
    mockery::stub(copy_file_to_bucket, "file.exists", function(x) TRUE)
    mockery::stub(copy_file_to_bucket, "Sys.which", function(x) "")
    expect_error(
        copy_file_to_bucket("local.txt", "gs://bucket", "out.txt"),
        "gsutil is not installed"
    )
})

test_that("copy_file_to_bucket calls gsutil cp and ls on success", {
    calls <- character(0)
    mockery::stub(copy_file_to_bucket, "file.exists", function(x) TRUE)
    mockery::stub(copy_file_to_bucket, "Sys.which",
                  function(x) "/usr/bin/gsutil")
    mockery::stub(copy_file_to_bucket, "system", function(cmd, intern) {
        calls <<- c(calls, cmd)
        character(0)
    })
    expect_invisible(
        copy_file_to_bucket("local.txt", "gs://bucket", "out.txt")
    )
    expect_equal(length(calls), 2L)
    expect_true(grepl("gsutil cp", calls[1]))
    expect_true(grepl("gsutil ls", calls[2]))
})

test_that("copy_file_to_bucket errors when gsutil cp exits non-zero", {
    mockery::stub(copy_file_to_bucket, "file.exists", function(x) TRUE)
    mockery::stub(copy_file_to_bucket, "Sys.which",
                  function(x) "/usr/bin/gsutil")
    mockery::stub(copy_file_to_bucket, "system", function(cmd, intern) {
        out <- "CommandException: bucket not found"
        attr(out, "status") <- 1L
        out
    })
    expect_error(
        copy_file_to_bucket("local.txt", "gs://bucket", "out.txt"),
        "gsutil cp failed"
    )
})

test_that("copy_file_to_bucket errors when gsutil ls exits non-zero", {
    call_count <- 0L
    mockery::stub(copy_file_to_bucket, "file.exists", function(x) TRUE)
    mockery::stub(copy_file_to_bucket, "Sys.which",
                  function(x) "/usr/bin/gsutil")
    mockery::stub(copy_file_to_bucket, "system", function(cmd, intern) {
        call_count <<- call_count + 1L
        if (call_count == 1L) {
            # gsutil cp succeeds
            return(character(0))
        }
        # gsutil ls fails
        out <- "CommandException: No URLs matched"
        attr(out, "status") <- 1L
        out
    })
    expect_error(
        copy_file_to_bucket("local.txt", "gs://bucket", "out.txt"),
        "gsutil ls failed"
    )
})
