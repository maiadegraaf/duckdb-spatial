# This file is included by DuckDB's build system. It specifies which extension to load

# Extension from this repo

# Disable tests on MinGW for 1.1.3, will be fixed in 1.1.4
if (MINGW)
    set(DO_TESTS "")
else ()
    set(DO_TESTS "LOAD_TESTS")
endif()

duckdb_extension_load(spatial
        SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}
        INCLUDE_DIR ${CMAKE_CURRENT_LIST_DIR}/src/spatial
        ${DO_TESTS}
)
