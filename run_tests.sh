#!/bin/sh
set -eu

OPENCV_PREFIX="${OPENCV_PREFIX:-third_party/opencv-install}"
OPENCV_INCLUDE="$OPENCV_PREFIX/include/opencv4"
OPENCV_LIB="$OPENCV_PREFIX/lib"

if [ ! -d "$OPENCV_INCLUDE" ] || [ ! -d "$OPENCV_LIB" ]; then
    echo "OpenCV local install not found at $OPENCV_PREFIX" >&2
    echo "Run ./scripts/build_opencv.sh first, or set OPENCV_PREFIX to an OpenCV install prefix." >&2
    exit 1
fi

g++ -std=c++17 -O2 -I"$OPENCV_INCLUDE" -Ithird_party \
    cpp/sfrmat5.cpp cpp/test_sfrmat5.cpp \
    -L"$OPENCV_LIB" -lopencv_core -Wl,-rpath,"$OPENCV_LIB" \
    -o /tmp/test_sfrmat5
/tmp/test_sfrmat5
