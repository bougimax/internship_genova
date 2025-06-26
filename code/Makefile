EXECUTABLE = ./build/cdt

INPUT_FILE = Input_file/dataset/good_meshes/60100.off

BUILD_CMD = cmake --build build --config Release
RUN_CMD = $(EXECUTABLE) -ovlpc $(INPUT_FILE)

.PHONY: all build run

all: build run

build:
	$(BUILD_CMD)

run:
	$(RUN_CMD)

