EVE 109 Section Week 8
================

### 1. Row and column names

Although we've dealt with row and column names in previous sections, we haven't discussed them explicitly. Once again, we'll use the `Arthritis` dataset:

``` r
library(vcd)
```

    ## Loading required package: grid

``` r
data(Arthritis)
head(Arthritis)
```

    ##   ID Treatment  Sex Age Improved
    ## 1 57   Treated Male  27     Some
    ## 2 46   Treated Male  29     None
    ## 3 77   Treated Male  30     None
    ## 4 17   Treated Male  32   Marked
    ## 5 36   Treated Male  46   Marked
    ## 6 23   Treated Male  58   Marked

 

Notice each column has a name: "ID", "Treatment","Sex","Age","Improved". We can ask for the columnames using the `colnames` command:

``` r
colnames(Arthritis)
```

    ## [1] "ID"        "Treatment" "Sex"       "Age"       "Improved"

 

If we do a calculation for the column, the column name does not count. For example, this does not return a category that says "Treatment"

``` r
table(Arthritis$Treatment)
```

    ## 
    ## Placebo Treated 
    ##      43      41

 

We can do a similar thing with row names. Look what happens when we ask for a summary of this dataset:

``` r
summary(Arthritis)
```

    ##        ID          Treatment      Sex          Age          Improved 
    ##  Min.   : 1.00   Placebo:43   Female:59   Min.   :23.00   None  :42  
    ##  1st Qu.:21.75   Treated:41   Male  :25   1st Qu.:46.00   Some  :14  
    ##  Median :42.50                            Median :57.00   Marked:28  
    ##  Mean   :42.50                            Mean   :53.36              
    ##  3rd Qu.:63.25                            3rd Qu.:63.00              
    ##  Max.   :84.00                            Max.   :74.00

 

It makes sense to summarize Treatment, Sex, Age, and Improved in this way, but it makes no sense at all to summarize ID, right? So instead, we can assign the IDs to the row names and remove them from the data frame:

``` r
rownames(Arthritis) <- Arthritis$ID
head(Arthritis)
```

    ##    ID Treatment  Sex Age Improved
    ## 57 57   Treated Male  27     Some
    ## 46 46   Treated Male  29     None
    ## 77 77   Treated Male  30     None
    ## 17 17   Treated Male  32   Marked
    ## 36 36   Treated Male  46   Marked
    ## 23 23   Treated Male  58   Marked

``` r
rownames(Arthritis)
```

    ##  [1] "57" "46" "77" "17" "36" "23" "75" "39" "33" "55" "30" "5"  "63" "83"
    ## [15] "66" "40" "6"  "7"  "72" "37" "82" "53" "79" "26" "28" "60" "22" "27"
    ## [29] "2"  "59" "62" "84" "64" "34" "58" "13" "61" "65" "11" "56" "43" "9" 
    ## [43] "14" "73" "74" "25" "18" "21" "52" "45" "41" "8"  "80" "12" "29" "50"
    ## [57] "38" "35" "51" "54" "76" "16" "69" "31" "20" "68" "81" "4"  "78" "70"
    ## [71] "49" "10" "47" "44" "24" "48" "19" "3"  "67" "32" "42" "15" "71" "1"

``` r
Arthritis <- Arthritis[,-1] #What does this do?!
head(Arthritis)
```

    ##    Treatment  Sex Age Improved
    ## 57   Treated Male  27     Some
    ## 46   Treated Male  29     None
    ## 77   Treated Male  30     None
    ## 17   Treated Male  32   Marked
    ## 36   Treated Male  46   Marked
    ## 23   Treated Male  58   Marked

Now if we summarize the dataset, it makes a little more sense:

``` r
summary(Arthritis)
```

    ##    Treatment      Sex          Age          Improved 
    ##  Placebo:43   Female:59   Min.   :23.00   None  :42  
    ##  Treated:41   Male  :25   1st Qu.:46.00   Some  :14  
    ##                           Median :57.00   Marked:28  
    ##                           Mean   :53.36              
    ##                           3rd Qu.:63.00              
    ##                           Max.   :74.00

 

Rownames and columnames can be useful to append identification tags to data when you don't want them to be part of any calculations. For example, when you have a matrix of genotypes you are working with, you can add sample names and marker IDs but maintain your ability to use apply functions across the entire dataframe. You will also find that some packages require row and column names in order to match dataframes.

     

### 1. GWAS

This week we'll use the data from the Bosse et al. (2017) paper to learn how to perform a Genome Wide Association Study (GWAS). There are many different softwares available for GWAS. We will use `rrBLUP` because it has an R implementation and a fairly straightforward file format. Install the package and load the library. Take a look at the help page for the command `GWAS`

``` r
install.packages("rrBLUP")
```

``` r
library(rrBLUP)
?GWAS
```

 

So, it looks like we need two data frames, one with phenotype data and one with genotype data. I have given you some example files. Notice how I use the arguement `row.names=1` to say that the first column should be used as the row names. This is because `rrBLUP` requires rownames for its calculations. Interestingly, R assumes you want column names unless you tell it otherwise.

``` r
pheno <- read.csv("data/example_pheno.csv",row.names=1)
head(pheno)
```

    ##              ID BillLength
    ## E509370 E509370      0.096
    ## F464559 F464559      0.065
    ## K578717 K578717      0.019
    ## L312169 L312169     -0.116
    ## L312889 L312889      0.221
    ## L315203 L315203      0.296

``` r
geno <- read.csv("data/example_geno.csv",row.names=1)
head(geno[,1:10])    # Just look at the first 10 columns
```

    ##                           snp chr   pos E509370 F464559 K578717 L312169
    ## AX.100678681_C AX.100678681_C   1   511       0       0       0      -1
    ## AX.100747228_A AX.100747228_A   1  5379      -1      -1      -1       0
    ## AX.100836047_C AX.100836047_C   1 21413      -1       0       0       0
    ## AX.100918667_G AX.100918667_G   1 23383      -1       0       0       0
    ## AX.100974310_A AX.100974310_A   1 26112       0       1       0      -1
    ## AX.100811471_T AX.100811471_T   1 28471       0       1       1      -1
    ##                L312889 L315203 L315204
    ## AX.100678681_C      -1      -1      -1
    ## AX.100747228_A      -1       1      -1
    ## AX.100836047_C      -1       0      -1
    ## AX.100918667_G       0       1       0
    ## AX.100974310_A       0       0       0
    ## AX.100811471_T       0      -1      -1

 

Notice that the genotype file format is a little different than what we've used before. You can see the description in on the `GWAS` help page. The first column (and the row names) are the SNP name, the second column is the chromosome, and the third column is the position. Each column after that is an individual and genotypes are coded as -1 (homozygote), 0 (heterozygote), or 1 (homozygote). Now we can use the genotype and phenotype files to perform the GWAS.

``` r
g <- GWAS(pheno,geno,min.MAF=0.05,plot=TRUE)
```

    ## [1] "GWAS for trait: BillLength"
    ## [1] "Variance components estimated. Testing markers."

 

How do we assess significance of these SNPs? Usually, we want to show a p-value cutoff. p&lt;0.05 is commonly used, as is p&lt;0.01. Notice that the y-axis on this plot is -log(p). We can add lines showing different p-value cutoffs like this:

``` r
g <- GWAS(pheno,geno,min.MAF=0.05,plot=TRUE)
```

    ## [1] "GWAS for trait: BillLength"
    ## [1] "Variance components estimated. Testing markers."

``` r
abline(h=-log(0.05),col="red",lty="dashed")
```

 

In this case none of our markers are significant. Let's take a look at the results generated by the `GWAS` command:

``` r
head(g)
```

    ##                           snp chr   pos BillLength
    ## AX.100678681_C AX.100678681_C   1   511 0.22120646
    ## AX.100747228_A AX.100747228_A   1  5379 0.05446761
    ## AX.100836047_C AX.100836047_C   1 21413 0.47992555
    ## AX.100918667_G AX.100918667_G   1 23383 0.02321833
    ## AX.100974310_A AX.100974310_A   1 26112 0.55926248
    ## AX.100811471_T AX.100811471_T   1 28471 0.01055644

 

Here the first column is the SNP ID, the second column is the chromosome, the third column is the position along that chromosome, and the fourth column is the association between bill length and genetic variation at that position (-log(p)).

     

Homework
========

Now that we know how to perform a simple GWAS, we can use the full data set from Bosse et al. (2017).

### *Homework 8: Write a script that does the following:*

#### 1. Load in the genotype data (snp\_data.csv) and the phenotype data (bill\_data.csv). Notice that these data are not in the same order.

#### 2. Perform and plot the GWAS for bill length. This may take a while. Add a dashed line representing a p-value cutoff of p&lt;0.05

#### 3. Using the GWAS output, make a barplot showing how many *significant* (p&lt;0.05) SNPs are on each chromosome (this will require application several skills we've learned previously)
