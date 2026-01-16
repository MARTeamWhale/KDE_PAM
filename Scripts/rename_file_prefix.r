# rename_file_prefix.R
# Removes leading "tif" from filenames in a folder, across all extensions.
# Example: tifKDE_foo.shx -> KDE_foo.shx

# ---- SET THIS ----
folder <- "/Users/chirp/CODE/KDE_PAM/output/shapes/baseline_match_v3/multi"   # change me

# ---- OPTIONS ----
prefix <- "tif"
dry_run <- F   # set to FALSE to actually rename

# ---- RUN ----
stopifnot(dir.exists(folder))

files_full <- list.files(folder, full.names = TRUE, all.files = FALSE, no.. = TRUE)
files_base <- basename(files_full)

idx <- startsWith(files_base, prefix)
to_rename_full <- files_full[idx]
to_rename_base <- files_base[idx]

if (length(to_rename_full) == 0) {
  message("No files starting with '", prefix, "' found in: ", folder)
  quit(save = "no")
}

new_base <- sub(paste0("^", prefix), "", to_rename_base)
new_full <- file.path(folder, new_base)

# Preview table
preview <- data.frame(
  from = to_rename_base,
  to   = new_base,
  stringsAsFactors = FALSE
)
print(preview, row.names = FALSE)

# Collision check (target already exists)
collide <- file.exists(new_full)
if (any(collide)) {
  message("\nThese target filenames already exist. Fix or move them first; no renames performed for collisions:\n")
  print(preview[collide, ], row.names = FALSE)
}

# Only rename non-colliding files
ok <- !collide
if (!any(ok)) {
  message("\nAll proposed renames collide. Nothing to do.")
  quit(save = "no")
}

if (dry_run) {
  message("\nDry run only. Set dry_run <- FALSE to apply ", sum(ok), " renames.")
} else {
  success <- file.rename(to_rename_full[ok], new_full[ok])
  if (!all(success)) {
    message("\nSome renames failed:")
    print(preview[ok, ][!success, ], row.names = FALSE)
  } else {
    message("\nRenamed ", sum(success), " files.")
  }
}

