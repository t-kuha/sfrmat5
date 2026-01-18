#!/bin/sh
set -eu

g++ -std=c++17 -O2 -Ithird_party/eigen-3.4.0 cpp/sfrmat5.cpp cpp/test_sfrmat5.cpp -o /tmp/test_sfrmat5
/tmp/test_sfrmat5
