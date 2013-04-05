ldpc:
	@echo del ldpc.jar if exist
	rm ldpc.jar -f
	@echo build 
	javac -s src -d build src/com/jove/ldpc/*.java
	@echo jar
	jar cf ldpc.jar -C build .


