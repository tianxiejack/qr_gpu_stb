################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/FindFtp.cpp \
../src/MotionFilter.cpp \
../src/kalman.cpp \
../src/matchingPoint.cpp \
../src/motionCompensate.cpp \
../src/preprocess.cpp \
../src/stable.cpp 

CU_SRCS += \
../src/cuda.cu 

CU_DEPS += \
./src/cuda.d 

OBJS += \
./src/FindFtp.o \
./src/MotionFilter.o \
./src/cuda.o \
./src/kalman.o \
./src/matchingPoint.o \
./src/motionCompensate.o \
./src/preprocess.o \
./src/stable.o 

CPP_DEPS += \
./src/FindFtp.d \
./src/MotionFilter.d \
./src/kalman.d \
./src/matchingPoint.d \
./src/motionCompensate.d \
./src/preprocess.d \
./src/stable.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/usr/include/opencv -I/usr/lib/gcc-cross/aarch64-linux-gnu/5/include -I../src/OSA_CAP/inc -I/usr/include -I/usr/include/opencv2 -I../inc -G -g -O0 -Xcompiler -fPIC -Xcompiler -fopenmp -ccbin aarch64-linux-gnu-g++ -gencode arch=compute_50,code=sm_50 -m64 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/usr/include/opencv -I/usr/lib/gcc-cross/aarch64-linux-gnu/5/include -I../src/OSA_CAP/inc -I/usr/include -I/usr/include/opencv2 -I../inc -G -g -O0 -Xcompiler -fPIC -Xcompiler -fopenmp --compile -m64 -ccbin aarch64-linux-gnu-g++  -x c++ -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/%.o: ../src/%.cu
	@echo 'Building file: $<'
	@echo 'Invoking: NVCC Compiler'
	/usr/local/cuda-8.0/bin/nvcc -I/usr/include/opencv -I/usr/lib/gcc-cross/aarch64-linux-gnu/5/include -I../src/OSA_CAP/inc -I/usr/include -I/usr/include/opencv2 -I../inc -G -g -O0 -Xcompiler -fPIC -Xcompiler -fopenmp -ccbin aarch64-linux-gnu-g++ -gencode arch=compute_50,code=sm_50 -m64 -odir "src" -M -o "$(@:%.o=%.d)" "$<"
	/usr/local/cuda-8.0/bin/nvcc -I/usr/include/opencv -I/usr/lib/gcc-cross/aarch64-linux-gnu/5/include -I../src/OSA_CAP/inc -I/usr/include -I/usr/include/opencv2 -I../inc -G -g -O0 -Xcompiler -fPIC -Xcompiler -fopenmp --compile --relocatable-device-code=false -gencode arch=compute_50,code=compute_50 -gencode arch=compute_50,code=sm_50 -m64 -ccbin aarch64-linux-gnu-g++  -x cu -o  "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


