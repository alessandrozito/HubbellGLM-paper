# Download the annual Human Footprint rasters from figshare.
#
# Produces: cache/hfp/

source(file.path(path.expand("~/HubbellGLM-paper"),
                 "preprocessing", "paths.R"))

if (!require("rfigshare")) {
  install.packages(
    "https://cran.r-project.org/src/contrib/Archive/rfigshare/rfigshare_0.3.8.tar.gz",
    repos = NULL, type = "source"
  )
  # alternative: remotes::install_github("ropensci/rfigshare")
  library(rfigshare)
}


# Accessing Data ----------------------------------------------------------

# The Human Footprint collection. It holds one zip per year (hfp2000.zip ...
# hfp2018.zip), each containing a single GeoTIFF.
article_id <- "16571064"

# Get the details for all the files contained in the article ID
metadata <- fs_details(article_id, mine = FALSE)

# Get the urls for download of each of the files
file_urls <- fs_download(article_id)

# Get the names of each of the files
file_names <- sapply(metadata$files, '[[', 'name')

# Combine URLs with files names to create a data frame
files_data <- data.frame(name = file_names,
                         url  = file_urls,
                         stringsAsFactors = FALSE)


# Subsetting --------------------------------------------------------------

# Subset to the year(s) you need. The sample runs 2007-2026 but the product
# stops at 2024, so later years carry no human footprint at all.
pattern <- "hfp2014"
files_data_subset <- files_data[grep(pattern, files_data$name), ]

# patterns <- paste0("hfp", 2010:2016)
# files_data_subset <- files_data[grep(paste(patterns, collapse = "|"), files_data$name), ]

stopifnot(nrow(files_data_subset) > 0)


# Downloading -------------------------------------------------------------

dir_name <- file.path(CACHE, "hfp")
dir.create(dir_name, recursive = TRUE, showWarnings = FALSE)

# Figshare throttles anonymous traffic and answers with a 403 whose body is a
# 118-byte HTML page. download.file() writes that to disk without complaining,
# so the .zip looks fine until unzip() fails with a confusing error. Check the
# PK\003\004 magic bytes before trusting the file.
is_valid_zip <- function(path) {
  if (!file.exists(path) || file.size(path) < 1024) return(FALSE)
  identical(readBin(path, "raw", n = 4), as.raw(c(0x50, 0x4b, 0x03, 0x04)))
}

options(timeout = 3600)  # the 2014 file is ~436 MB

invisible(lapply(seq_len(nrow(files_data_subset)), FUN = function(i) {
  zip_path <- file.path(dir_name, files_data_subset[i, "name"])

  download.file(url      = files_data_subset[i, "url"],
                destfile = zip_path,
                mode     = "wb")

  if (!is_valid_zip(zip_path)) {
    unlink(zip_path)
    stop(files_data_subset[i, "name"], " is not a zip archive - figshare ",
         "returned an error page instead of the data. See the note below.")
  }

  unzip(zipfile = zip_path, exdir = dir_name)
}))

list.files(dir_name)


# If you get 403 Forbidden ------------------------------------------------

# A 403 from figshare is not an R problem: every route (rfigshare,
# download.file, httr, curl, or a browser) requests the same figshare hosts,
# and figshare's load balancer can reject a whole IP range before any
# authentication is evaluated. Two ways to tell them apart:
#
#   curl -sI https://figshare.com/
#
# If plain figshare.com also returns 403 with "server: awselb/2.0", the
# network itself is blocked and no token or package will help - retry from a
# different network (home vs. institutional, or a VPN), or download the zip in
# a browser and drop it in data/cache/hfp/.
#
# If figshare.com loads but the download 403s, it is rate limiting. Create a
# free personal token at figshare.com > Account settings > Applications >
# Create personal token, add FIGSHARE_PAT=<token> to ~/.Renviron, restart R,
# and use httr directly:
#
#   httr::GET(url, httr::add_headers(
#     Authorization = paste("token", Sys.getenv("FIGSHARE_PAT"))),
#     httr::write_disk(zip_path, overwrite = TRUE), httr::progress())
