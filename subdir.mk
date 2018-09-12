################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../EquityTests.cpp \
../HUNLtests.cpp \
../IntegrationTests.cpp \
../stdafx.cpp 

OBJS += \
./EquityTests.o \
./HUNLtests.o \
./IntegrationTests.o \
./stdafx.o 

CPP_DEPS += \
./EquityTests.d \
./HUNLtests.d \
./IntegrationTests.d \
./stdafx.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Intel C++ Compiler'
	icpc -O3 -qopt-prefetch=2 -inline-level=2 -I"/home/brans/repos/deepstack/DeepStackCpp" -I/usr/local/cuda/include/ -I"/home/brans/repos/deepstack/IntegrationTests/../OMPEval" -I"/home/brans/repos/deepstack/IntegrationTests/../eigen" -I"/opt/intel/compilers_and_libraries_2018/linux/mkl/include" -I"/home/brans/repos/deepstack/IntegrationTests/../DeepStackCpp/" -I"/home/brans/repos/deepstack/IntegrationTests/../ziplib/Source/ZipLib" -DUSE_GPU -std=c++14 -fp-model fast=2 -mtune=skylake -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -c -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


