1. What is the structure of each of the three input files for SiSaNA's preprocessing step? What do the columns represent in each of the files? (Hint: you can use the 'head' command in bash to view just the first few lines of a file, and find additional information on the file inputs inside the params.yml file)

2. After running the generate step, a few files are automatically generated. Which of them might you use to determine which genes were being most highly regulated?

3. Following the comparison of the expression of the two Luminal breast cancer groups (the default setup in the params.yml file), take a look at the output .txt file. What are some of the most differentially expressed genes? 

4. Do the same comparison as before but for the indegrees this time. 

5. Take a look at the output of the "sisana visualize survival" command. Do you know how to interpret this image? Additionally, there is a p-value listed there. Do you think we can conclude based on this p-value that there is a significant difference in the survival of Luminal A vs Luminal B patients?

6. Following running the "sisana compare gsea" step, what are the three top pathways? Run this step again, but with the Hallmark gene set instead of the Reactome one.

7. Now let's look at the output of the "sisana visualize volcano" command. What do the dots on the left and right side of the plot represent? Do you recognize any of these genes?

8. Take 5 of the top differential indegrees (which you calculated in step 4) and make a box plot for those genes.

9. How many TFs have edges that connect to (regulate) the GDI2 gene in your network? (Hint: Use the help function of SiSaNA (sisana -h) to find a subcommand that can help with this).
