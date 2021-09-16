# dzd_assessment
Coding assessment for Day Zero Diagnostics

### "How did you test the code for each problem?"

For problem 1a (count all the 21-mers within the sequences), I tested in two ways. First, I wanted to ensure that my indices were correct and that I was considering each kmer in the sequence (not leaving any off). To do this, I simply used small test cases and manually confirmed that I was seeing all the kmers I expected. Then, to ensure that my frequency counts were accurate, I added a few lines to the function problem1a() which computed the expected number of frequencies and the actual number.

For problem 1b (working with SQLite), I simply utilized the testing commands described in the provided SQLite tutorial. I included a few of those lines of code, commented out, in the function problem1b(). Trusting that my initial k-mer frequency table was accurate (I had no reason to believe it wasn't based on the testing from above), I further trusted that the SQLite command accurately added the information into the kmer table (...something something Law of Commutativity :D).

Problem 2 was tricky to test, as it's tedious and unrewarding to manually check for all possible kmers in a sequence and see if they are within 2 letters' difference. So, as I did in testing problem 1a, I tried a few small examples (including the one provided in the assessment) and verified that I saw what I expected. In terms of testing the running time, I wrote a few lines of code in function problem2() which generate random sequences of arbitrary length and runs the code on them.

### "How will the code perform with a very large input file?"

Problem 1 runs in linear time - O(n). We first read the file into memory (one pass through the file), and then we read through each sequence once (for a combined second pass through the file). So, as the file gets larger, the time will get slower in a directly proportional manner. If the file gets so large that we can't load it into memory, then the code will not work.

Problem 2 runs in O(n^2) time. As we run the code, one pass through the sequence is performed because we consider each k-mer in seq. However, we are also comparing each k-mer to kseq, which has a size with an upper bound equal to the length of seq. So, we are essentially running in a time equal to O(k * n), where k is the length of kseq and n is the length of seq. But, since k could equal n, this is the same as O(n^2). So, as the file get larger, the running time will increase in an exponential manner.

### How to run the code:

main() contains three function calls, one for each of problem 1a, problem 1b, and problem 2. You can modify which of these you would like to run by commenting/uncommenting these three lines. Note that there is no need to run the function for problem 1b, as the table is already created. Simply issue the command-line commands included in that function to see the table for yourself. All this to say, running "python3 solutions.py" once you have modified main() should do the trick :)
