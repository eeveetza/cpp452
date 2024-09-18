#!/bin/bash
 
outdir=build
 
## cleanup
 
# rm -rfd $outdir
# rm -rfd bin
 
## build
 
buildType=$1
 
if [[ $buildType != Debug ]]; then
    buildType=Release
fi
 
enableCodeCoverage=OFF
 
if [[ $2 == Coverage ]]; then
    enableCodeCoverage=ON
fi
 
echo "cmake -DCMAKE_BUILD_TYPE=$buildType -S . -B $outdir"
cmake -DCMAKE_BUILD_TYPE=$buildType -DENABLE_CODE_COVERAGE:BOOL=$enableCodeCoverage -S . -B $outdir
cmake --build $outdir 2>outputdebug.txt
grep -n "error:" outputdebug.txt
  
## test
 
ctest --test-dir $outdir
 
# pushd bin
# ./$executable
# popd
 
# if [[ $2 == Coverage ]]; then
#     lcov -c -d . -o code-coverage.info
#     genhtml code-coverage.info --output-directory coverage-report-html
# fi