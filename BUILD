load("@rules_foreign_cc//foreign_cc:defs.bzl", "cmake")

filegroup(
    name = "srcs",
    srcs = glob(["**"]),
    visibility = ["//visibility:public"],
)

cmake(
    name = "gridlabd",
    build_args = [
        "-j 8",
    ],
    lib_source = ":srcs",
    out_binaries = ["gridlabd"],
    out_data_dirs = ["share"],
    targets = ["install"],
)
