EXECUTABLE = ./build/cdt

INPUT_FILE = Input_file/112544.off

BUILD_CMD = cmake --build build --config Release
RUN_CMD = $(EXECUTABLE) -o $(INPUT_FILE)

.PHONY: all build run

all: build run

build:
	$(BUILD_CMD)

run:
	$(RUN_CMD)

