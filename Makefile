src:=$(wildcard src/com/jove/ldpc/*.java) 
ldpc.jar:${src}
	@echo build 
	javac -s src -d build src/com/jove/ldpc/*.java
	@echo jar
	jar cf ldpc.jar -C build .

build/LDPCTest.class:ldpc.jar test/LDPCTest.java
	javac -cp ./ldpc.jar -d build test/*.java
	
test:build/LDPCTest.class
	java -cp build:ladpc.jar LDPCTest examples/DenEvl_9.deg examples/DenEvl_9.code 0.25 28 72 83

