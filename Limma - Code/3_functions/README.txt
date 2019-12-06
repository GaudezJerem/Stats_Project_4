var_names

Description
-------------

Extract substrings in a character vector that have varying positions in each string element in the vector

Usage
-------------

colData, text
		a character vector of sample names

pattern
		regular expression 
		Provide pattern within "x"

begin
		integer. Position relative to pattern of first chr of substr to be extracted
		boolean. Set to FALSE if first chr is at the start of the string element for all substr to be extracted

end
		integer. Position relative to pattern of last chr of substr to be extracted

Details
-------------

The function uses the function regexpr to return a vector of integers corresponding to the first and last chr of the
substring to be extracted, as dictated by the begin and end arguements relative to the pattern given.

The function then uses the function substr to extract the desired substring via the vector of integers dictating 
the start and stop argument values for this function.

By setting the begin arguement to FALSE the substr start value will be set to 1 and the first chr of each str 
element will be returned.

Value
-------------

Character vector of extracted substrings.

Examples
-------------
1. With boolean

colData_example <- c("Mock_t09_rep1", "Infected_t09_rep1")

example <- var_names(colData = colData_example,
	             pattern = "_",
	             begin = FALSE,
	             end = -1)

head(example)
> "Mock", "Infected"

2. Without boolean

colData_example <- c("Mock_t09_rep1", "Infected_t12_rep1")

example <- var_names(colData = colData_example,
	             pattern = "_r",
	             begin = -2,
	             end = -1)
head(example)
> "09", "12"