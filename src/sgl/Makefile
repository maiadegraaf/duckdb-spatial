.PHONY: all clean report coverage

LLVM_PROFDATA = $(shell brew --prefix llvm)/bin/llvm-profdata
LLVM_COV = $(shell brew --prefix llvm)/bin/llvm-cov

all: sgl_test

sgl_test: sgl.hpp sgl.cpp sgl_test.cpp
	$(CXX) -o $@ sgl.cpp sgl_test.cpp -std=c++11 -g -O0 -fprofile-instr-generate -fcoverage-mapping -DSGL_COVERAGE -fsanitize=address -fsanitize=undefined

sgl_test.profdata: sgl_test
	LLVM_PROFILE_FILE="sgl_test.profraw" ./sgl_test
	$(LLVM_PROFDATA) merge -output=sgl_test.profdata sgl_test.profraw

# Print coverage report
report: sgl_test.profdata
	$(LLVM_COV) report sgl_test -instr-profile=sgl_test.profdata -use-color -show-region-summary=false

# Find lines with missing code coverage
coverage: sgl_test.profdata
	$(LLVM_COV) show sgl_test -instr-profile=sgl_test.profdata \
		-use-color \
		-show-region-summary=false \
		-show-line-counts=false \
		-line-coverage-lt=100

clean:
	rm sgl_test
	rm sgl_test.profraw
	rm sgl_test.profdata