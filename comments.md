
##  81 thoughts on “How to use DESeq2 to analyse RNAseq data”

  1. ![](./DESeq2_files/71ae912f875c12050f8af47f9bf3775e)Rupesh on [February 19, 2014 at 11:49 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-412) said:

Dear sir,  
Thanks sir for posting this great tutorial. Sir i would like to plot MA with
padj =< 0.05. I know its by default deseq2 set it as 0.1….

could you help me to plot RED dots (significant DE at padj =< 0.05) in plotMA.
I assumed i have to adjust this value in dds ???

Thanks and looking forward for your kind help

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=412#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [February 19, 2014 at 7:52 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-415) said:

Hi There,

PlotMA is a function from DESeq, you can see all the info here at the
biocondutor page listing all the doc material [http://bioconductor.org/package
s/2.13/bioc/manuals/DESeq/man/DESeq.pdf](https://web.archive.org/web/201603282
34446/http://bioconductor.org/packages/2.13/bioc/manuals/DESeq/man/DESeq.pdf).
I haven’t checked the DESeq2 function, but I guess it would be the same.

the function is listed as (from the above):  
plotMA(x, ylim,  
col = ifelse(x$padj>=0.1, “gray32”, “red3”),  
linecol = “#ff000080”,  
xlab = “mean of normalized counts”, ylab = expression(log[2]~fold~change),  
log = “x”, cex=0.45, …)

So, following that logic we just need to pass a col variable to override the
default, so this should work for your needs:

plotMA(dds,col=ifelse(res$padj>0.05,”grey32″,”red”))

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=415#respond)

      - ![](./DESeq2_files/71ae912f875c12050f8af47f9bf3775e\(1\))Rupesh on [February 20, 2014 at 10:33 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-419) said:

Thanks sir but logic is different….from deseq to deseq2  
i got it, it should be ….  
plotMA(dds,pvalCutoff=.05,ylim=c(-2,2),main=”DESeq2″)

      - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [February 20, 2014 at 7:02 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-421) said:

Thanks for posting that, good to know.

  2. ![](./DESeq2_files/7e5062602b15fe72aff056e9ea0f61a2)Will on [February 27, 2014 at 10:48 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-449) said:

Hi,  
I am having trouble getting started and it is probably due to the fact my
counts tables did not derive from HT-Seq (I used a different software, CLC). I
do, however, have a merged file containing the gene ID, followed by the counts
of each of the different samples. I started following your DESeq tutorial and
got to the point of  
conds<-factor(c(""untreated","treated"…..)) and then ran into problems because
DESeq2 doesn't have newCountDataSet function. This, of course, led me to this
tutorial but don't know how to proceed since I already created a merged count
file.

Apologies for the basic question. I'd appreciate any assistance you can
provide.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=449#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [February 27, 2014 at 11:58 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-450) said:

The simplest (but rather unsatisfying) solution is to just split the table
using a script or excel, col1 gene names, cold 2 counts for each replicate.

Going by this post ([http://seqanswers.com/forums/showthread.php?t=39253](http
s://web.archive.org/web/20160328234446/http://seqanswers.com/forums/showthread
.php?t=39253))though it looks like the DESeqDataSetFromMatrix should get you
started (“?DESeqDataSetFromMatrix”) to see usage

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=450#respond)

      - ![](./DESeq2_files/7e5062602b15fe72aff056e9ea0f61a2\(1\))Will on [February 28, 2014 at 2:12 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-451) said:

Thanks! I figured out how to get it in. Here’s the command lines:  
>countsTablerownames(countsTable)countsTablehead(countsTable)  
>countDatacolData<-data.frame(condition=factor(c("insert your treatments in
proper order for each column of countsTable")))  
ddscolData(dds)$condition<-factor(colData(dds)$condition,levels=c("your
treatments"))  
dds

This got me to the point where I can run DESeq and convert the dds file to
res. Now I am trying to figure out what is wrong with the PlotMA function. I
type:  
plotMA(dds,ylim=c(-2,2),main="DESeq2")

and receive the following error: 'x' must be a data frame with columns named
'baseMean', 'log2FoldChange'.

Not sure why I am getting this. Will need to troubleshoot.

      - ![](./DESeq2_files/7e5062602b15fe72aff056e9ea0f61a2\(1\))Will on [March 4, 2014 at 5:08 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-490) said:

Hi,

Just a general note for anyone going through this: make sure you have the
updated packages before beginning. I know I ran into some problems with this.

Thanks Dave for the great tutorial. I found it to be very helpful!

A couple of things: I had trouble with plotting the “h1 <\-
hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)". My R didn't know the
function of "use" or "!use". I compared sessionInfo() and didn't see anything
stand out. I'm sure I am overlooking something simple here…

Lastly, your line, resClean write.csv(as.data.frame(resClean),file="sim_condit
ion_treated_results_cleaned_deseq2.csv")

gave me the error

Error: unexpected symbol in "resClean write.csv"

Other than that, very straightforward. Thanks again!

      - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 5, 2014 at 5:55 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-496) said:

Thanks for the thanks.

I’ve been having trouble with the code being mushed up by wordpress. I was
making a stupid mistake here, that hopefully is fixed. I’ve re-entered all the
code (probably adding in some mistakes along the way) and the formating looks
a lot better. Hopefully this will fix some of these silly errrors like
“”resClean write.csv” which is just complete mushed code.

  3. ![](./DESeq2_files/fcc61222eec44ced29caa22fea7df5d1)Super on [March 3, 2014 at 12:00 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-475) said:

Hi Dave，nice post.  
Is edgeR/DESeq/DESeq2 only outputs the list of genes? If I want to find the DE
transcripts or coding sequence CDS, could I do it by DESeq?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=475#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 3, 2014 at 7:34 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-479) said:

Thanks for the comment! If you have a well annotated GTF file that contains
isoform information and a splice aware mapper like tophat HTSeq count will use
that information to generate isoform specific read counts (see [http://www-hub
er.embl.de/users/anders/HTSeq/doc/count.html](https://web.archive.org/web/2016
0328234446/http://www-huber.embl.de/users/anders/HTSeq/doc/count.html)). My
bug of interest does not have good gene models at the moment so I’m stuck with
just using loci level information in all the example tables I put up here.
Other options are to try DEXseq which was designed to deal with alt splice
forms, but personally I have never used this. I think its correct to say the
biggest difference between cufflinks and DESeq is that the former will
(additionally) use the reads to generate novel gene models that are then
merged with the genome gene models.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=479#respond)

      - ![](./DESeq2_files/fcc61222eec44ced29caa22fea7df5d1\(1\))Super on [March 4, 2014 at 11:01 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-484) said:

Thank you! Could I summarize the Cufflinks vs DESeq  
1.Cufflinks is stronger for their efficiency and test rate (like the TPR?).?  
2.Cufflinks could potentially discovering novel genes and transcripts? and
output the DE transcripts/CDS/TSS which DESeq probably not (decided by HTseq
output).  
3\. DESeq is more suitable for unbalanced-samples (like 3 samples in treated
and 4 samples in untreated) because Cufflinks/Cuffdiff is only suitable for
pairwise comparison.  
Am I right? Do you have your revised/supplement comments?

      - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 5, 2014 at 6:04 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-497) said:

I think its fair to say that cufflinks2 (the cuffdiff part) and deseq are now
more similar than they used to be with regard to the DEG stats (see [http://cu
fflinks.cbcb.umd.edu/manual.html#dispersion_meth](https://web.archive.org/web/
20160328234446/http://cufflinks.cbcb.umd.edu/manual.html#dispersion_meth)).
But the other big difference is normalisations based on FPKM versus counts,
its perhaps an unsettled debate which is best.

      - ![](./DESeq2_files/fcc61222eec44ced29caa22fea7df5d1\(1\))Super on [March 5, 2014 at 12:06 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-503) said:

Thank you. I don’t know counts or FPKM which is better as well.So what do you
think? Compare the edgeR/DESeq/Cuffdiff? Do you have any suggestion or
conclusions? (e.g. running time? statistic power? and so on) Thank you!

      - ![](./DESeq2_files/7e5062602b15fe72aff056e9ea0f61a2\(1\))Will on [March 5, 2014 at 4:08 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-505) said:

I thought I would chime in my thoughts on this topic… There was a fairly
recent (2013) paper that compares Cufflinks2 and DESeq2 (and EdgeR and a
couple other methods). One of the conclusion of this study suggested Cufflinks
FPKM may be a bit too conservative and overlooks many differential gene
expressions. Personally, I think it is prudent to try out a couple of these
programs and draw conclusions from both. For example, if both a conservative
and liberal analysis produce the same set of genes as being differentially
expressed then you have greater evidence than if only one analysis produces
the results.

      - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 7, 2014 at 10:31 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-525) said:

I generally agree, although I always throught FPKM were actually too liberal
at least based in this (excellent) paper [http://bib.oxfordjournals.org/conten
t/early/2012/09/15/bib.bbs046.abstract](https://web.archive.org/web/2016032823
4446/http://bib.oxfordjournals.org/content/early/2012/09/15/bib.bbs046.abstrac
t). But yes trying different methods is important. I think we have to be
realistic and remember that we are asking a lot of the models to accurately
estimate expression based on just 3 samples. I guess when seq gets even
cheaper we can start using many more biological replicates to get a better
sample of the transcriptome under some condition or at least make sure we have
enough coverage to get good readings for low expressed genes that are at the
edge of the models statistical power.

  4. ![](./DESeq2_files/7e5062602b15fe72aff056e9ea0f61a2)Will on [March 6, 2014 at 12:07 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-508) said:

Hi again,  
I had a general question in interpreting the final results. I noticed in my
data I have a long list of genes with an ever-increasing p-value. In fact, I
don’t start seeing “NA” under “padj” until row 22685. Is this to say that I
have 22,684 significant differential expressed genes? If not, how can I
justify statistically where to draw the minimum p-value threshold? I was
looking through the DESeq2 manual, and although you don’t plot it here, I was
wondering the Benjamini-Hochberg multiple testing adjustment procedure (Fig.12
pg 30) could be used to reduce the number of genes to examine. If I understand
it correctly, anything to the left of the red and black line intersection is
controlled for FDR at alpha 0.01. Thus, only 1481 genes (in the manual data)
are of interest.

Your thoughts?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=508#respond)

  5. ![](./DESeq2_files/44fdfe69e9c82484fca3165982f5c7c7)geneart on [March 26, 2014 at 12:06 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-624) said:

Thanks for a great post Dr.Wheeler ! Very much appreciate your time and effort
in doing this which is a wonderful guide for a rookie like myself ! As I was
following your analysis steps , I got an error while plotting plotMA. I had
this same error while running DESeq as well.If I plot the dispersions against
counts it works but not by the way you do.THis is what it says:  
> plotMA(dseq,pvalCutoff=.05,ylim=c(-2,2))  
Error in as.vector(data) :  
no method for coercing this S4 class to a vector  
I looked it up on the net , and it seems like this is an R related errors and
not DEseq per sayt. However I am not sure. Please could you help me
troubleshoot this?  
THanks very much in advance !  
geneart.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=624#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 26, 2014 at 3:53 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-626) said:

Hi There, thanks for the kind words. Bit hard to work out whats wrong, but one
easy check is to make sure you are using the same versions as me by looking at
my sessionInfo() versus yours. Better yet, you can access the vignettes for
the version of DESeq you are using with this, browseVignettes(package =
“DESeq”) then double check that function in the doc to see if it has changed
at all. That would be a good place to start as often the functions change
slightly between versions.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=626#respond)

      - ![](./DESeq2_files/44fdfe69e9c82484fca3165982f5c7c7\(1\))geneart on [March 26, 2014 at 5:14 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-628) said:

Thanks for the reply Dr.Wheeler. I checked on my sessionInfo and nothing seems
to be drastically different from your sessionInfo, except that I am using a
different version of R , not too far off though ! also I checked the vignette
that shows the same function as u mention, so no problems there:  
> sessionInfo()  
R version 3.0.1 (2013-05-16)  
Platform: i386-w64-mingw32/i386 (32-bit)

locale:  
[1] LC_COLLATE=English_United States.1252  
[2] LC_CTYPE=English_United States.1252  
[3] LC_MONETARY=English_United States.1252  
[4] LC_NUMERIC=C  
[5] LC_TIME=English_United States.1252

attached base packages:  
[1] grDevices datasets parallel stats graphics utils methods  
[8] base

other attached packages:  
[1] DESeq2_1.2.10 RcppArmadillo_0.4.100.2.1  
[3] Rcpp_0.11.1 GenomicRanges_1.14.4  
[5] XVector_0.2.0 IRanges_1.20.7  
[7] DESeq_1.14.0 lattice_0.20-15  
[9] locfit_1.5-9.1 Biobase_2.22.0  
[11] BiocGenerics_0.8.0 BiocInstaller_1.12.0

loaded via a namespace (and not attached):  
[1] annotate_1.40.1 AnnotationDbi_1.24.0 DBI_0.2-7  
[4] genefilter_1.44.0 geneplotter_1.40.0 grid_3.0.1  
[7] RColorBrewer_1.0-5 RSQLite_0.11.4 splines_3.0.1  
[10] stats4_3.0.1 survival_2.37-7 tools_3.0.1  
[13] XML_3.98-1.1 xtable_1.7-3

However I did look for plotMA problems on web and one suggestion I got from
the DESeq2 manual is to detach BiocGenerics or load limma package later as
they both use the plotMA and somehow it interferes. I tried detaching the
BiocGenerics but this is what I got:  
> detach(“package:BiocGenerics”, unload=TRUE)  
Error: package ‘BiocGenerics’ is required by ‘GenomicRanges’ so will not be
detached  
> detach(“package:limma”, unload=TRUE)  
Error: package ‘limma’ is required by ‘edgeR’ so will not be detached  
> detach(“package:edgeR”, unload=TRUE)  
Warning messages:  
1: In FUN(X[[2L]], …) :  
Created a package name, ‘2014-03-26 12:55:56’, when none found  
2: In FUN(X[[2L]], …) :  
Created a package name, ‘2014-03-26 12:55:56’, when none found  
3: In FUN(X[[2L]], …) :  
Created a package name, ‘2014-03-26 12:55:56’, when none found  
4: In FUN(X[[2L]], …) :  
Created a package name, ‘2014-03-26 12:55:56’, when none found  
> detach(“package:limma”, unload=TRUE)  
There were 18 warnings (use warnings() to see them)  
> plotMA()  
Error in class(object) : ‘object’ is missing  
> plotMA(dseq,pvalCutoff=.05,ylim=c(-2,2))  
Error in (function (classes, fdef, mtable) :  
unable to find an inherited method for function ‘extractROWS’ for signature
‘”list”’

So at this point I am lost ! I would love to have this plotMA working for me ,
but cannot figure out the problem. Any suggestions? anyone ?Thanks :)  
geneart.

      - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 26, 2014 at 7:29 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-630) said:

The order is important, the last package to load will overwrite the previous
function. I wouldn;t detact the other other library and DESeq needs it. I
would load both and then explicitly call the DESeq2 function like this.

DESeq2::plotMA(bla bla)

#replacing bla bla with arguments.

Sorry I not really an R expert, but hopefully that helps

  6. ![](./DESeq2_files/44fdfe69e9c82484fca3165982f5c7c7)geneart on [March 26, 2014 at 10:02 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-631) said:

Thanks , will try this ang update again :)  
Geneart

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=631#respond)

  7. ![](./DESeq2_files/44fdfe69e9c82484fca3165982f5c7c7)geneart on [April 1, 2014 at 10:43 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-670) said:

Hello Dr.Wheeler,

> ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
directory=directory, design=~sampleCondition)  
Error in `row.names<-.data.frame`(`*tmp*`, value = value) :  
duplicate 'row.names' are not allowed  
In addition: Warning message:  
non-unique value when setting 'row.names': ‘balt.txt’

does this indicate I have a gene_name duplicated in the balt.txt file? If I
use HTSeqCOunt function instead of CountDataSet I get this error. If I use the
CountDataSet everything works fine. ANy suggestions? I am reading the files
just like how you did.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=670#respond)

  8. ![](./DESeq2_files/44fdfe69e9c82484fca3165982f5c7c7)geneart on [April 12, 2014 at 8:05 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-724) said:

Hello Dr.Wheeler,  
Finally after some digging through I got to figure how to fix my plotMA and it
is working. However I had one question though. In your demo you have used the
treated first and then you control in the sample conditions. However in
calling the level , later in ddsHTSEq u have called
levels=c(“untreated”,”treated”). I believe thi sshould nto be an issue.
However wanted to double check on this if it would be an issue??  
I have used control and treated as the order in my sample conditions an dhave
called levels=c(“control”,”treated”). when I plotted my
DESeq2::plotMA(dseq,pvalCutoff=.05,ylim=c(-2,2)) I see only 2 red spots.  
I was wondering if my order was an issue or just that I ahve two outliers ?
Please could you help me? HEre is my plot.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=724#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [April 16, 2014 at 10:42 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-751) said:

The order is only important for the direction of the fold change, it won’t
change any of the statistics. It makes sense to use control vs treatment, as
anything down in response to the treatment will be given a negative fold
change, and visa versa. Most people would expect a 2 fold up to equal UP in
response to treatment, for example. Only 2 significant DEG is better than
none!

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=751#respond)

  9. ![](./DESeq2_files/7d6ca16d4382bafce5f5825722bfc2bf.jpeg)Gthm on [April 14, 2014 at 3:40 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-727) said:

Should we filter the data before feeding in to DESEQ2 ? Like, remove entities
that have zero across all the conditions and then feed to DESEQ ?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=727#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [April 16, 2014 at 10:51 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-752) said:

Filtering the zero rows will certainly help reduce the multiple testing issues
by reducing the number of tests (by removing those that are guaranteed to
fail).  
Cheers

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=752#respond)

  10. ![](./DESeq2_files/e86645a0ba7ca9c6bee24953c60e199e)Ryan on [April 22, 2014 at 6:01 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-779) said:

Hi Dr. Wheeler-

I appreciate the time you’ve taken in producing this (and the previous DESeq)
tutorial. I had two questions regarding your python script, merge_tables.py,
for merging the count tables. (I employ HTSeq prior to running the
merge_tables.py.)

1) Prior to running your script (or before running DESeq), should I remove the
rows without gene counts? That is, the rows that contain no_feature, ambigous,
too_low_aQual, not_aligned, and alignment_not_unique?

2) Are there limitations on the number of count files that I can merge at one
time? And if so, would it be appropriate to run the script twice, on half the
count files and then merge them after the fact?

Cheers

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=779#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [April 22, 2014 at 9:54 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-781) said:

Hi Ryan,  
1) The script expects all the files to share the same row names (identifiers)
so don’t remove low/no count genes before you merge. You can do it after the
merge (even using excel if you don’t know scripting, just sort and chop the
rows that are zero across ALL columns). By zero counts I would mean zero
across all treatments, so you should only chop these once you see the results
for all the columns otherwise its not sensible thing to do.  
2) You can merge as many files as you like, just add the files and column
headings to the guide file. DONT run the script twice, it will break. Just add
all the files and corresponding column headings to the guide file and it will
work.  
Cheers  
Dave

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=781#respond)

  11. ![](./DESeq2_files/fcc61222eec44ced29caa22fea7df5d1)Quan on [May 19, 2014 at 10:02 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-933) said:

Hi Dave  
could you show me some scripts if I want to do multiple condition comparison?  
So far I have 3 conditions and 3 samples per condition, I don’t want to do
three times like C1 vs C2, C2 vs C3, C1 vs C3, could you please tell me how to
dow it just one time in DESeq2?  
Thank you!  
Quan

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=933#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [May 19, 2014 at 7:53 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-936) said:

I’m pretty sure you have to do them pairwise with DESeq. From memory Cufflinks
lets you do multiple comparisons at once, not sure.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=936#respond)

      - ![](./DESeq2_files/fcc61222eec44ced29caa22fea7df5d1\(1\))Quan on [May 20, 2014 at 10:54 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-943) said:

And from your case PCA plot(in this blog), which is not total separate two
classes(green and blue points), is this means that the difference expression
between two conditions is not obvious?  
Furthermore , is it similiar to the edgeR’s MDS plot?

      - ![](./DESeq2_files/fcc61222eec44ced29caa22fea7df5d1\(1\))Quan on [May 20, 2014 at 10:56 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-944) said:

PS: what is difference bwtween MDS and PCA plot meaning?

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [May 20, 2014 at 7:57 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-945) said:

Quan, correct, there is no grouping of treatments, but this was the result I
was hoping for because it was a comparison between two controls (see the text
above the figure). I;m not an expert like Google(-; the top hit for your other
question looks pretty good: [http://stats.stackexchange.com/questions/14002
/whats-the-difference-between-principal-components-analysis-and-multidimension
al](https://web.archive.org/web/20160328234446/http://stats.stackexchange.com/
questions/14002/whats-the-difference-between-principal-components-analysis-
and-multidimensional)

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=945#respond)

      - ![](./DESeq2_files/fcc61222eec44ced29caa22fea7df5d1\(1\))Quan on [May 21, 2014 at 8:51 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-946) said:

Cheers!

  12. ![](./DESeq2_files/24f771b1b82aae4d052739fc0593e5af)Julia on [May 25, 2014 at 6:33 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-971) said:

Hi!  
First of all, thank you for the great tutorial!  
However, I have a problem now with
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
directory=directory, design=~condition)

At the console,  
Error in file(file, "rt") : cannot open the connection  
In addition: Warning message:  
In file(file, "rt") :  
cannot open file 'C/Users/Julia/Desktop/I_375_Ago-IP.gene.hts': No such file
or directory  
appears.  
Can you help me?

Best wishes  
Julia

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=971#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [May 25, 2014 at 9:44 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-972) said:

Hi Julia, your file path looks wrong. Try:  
directory<-file.path('C:','Users/Julia/Desktop')

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=972#respond)

      - ![](./DESeq2_files/24f771b1b82aae4d052739fc0593e5af\(1\))Julia on [May 26, 2014 at 6:04 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-973) said:

thanks a lot! Now it works! :-)

  13. ![](./DESeq2_files/71ae912f875c12050f8af47f9bf3775e)Rupesh on [June 19, 2014 at 8:49 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1098) said:

Dear sir,

hello again

just few confusion need to resolve from your great experience.

Just wondering to compare DESeq2 and edgeR results of DE. But confusion is
that…..In edgeR it is advisable to filter low tags by cpm counts before DE
analysis while in DESeq2 its not. while The package DESeq2 consider
independent fileting (via function results ) step that is not available for
edgeR ??

so my query is that  
1) it is necessary to take filtration step of edgeR in edgeR (or its just
optional)???  
#####case 1 ( followed default functions)  
2) while comparing the results of DE gene of DESEq2 (without via edgeR of cpm
just used default functions i.e. deseq and results ) and in edgeR (with
filtration via cpm), i got results almost similar and that is what, i was
expecting….  
####case 2 (applied filtration in both via cpm )  
but while applying same filtering of edgeR via cpm in DESeq2 and then used the
function deseq and results (with and without independent filtering by both;
just for check), results are not similar et all…..very far DE in both (taking
0.05 padj as cutoff).

So please clear my doubt and suggest me what to do ???

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1098#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [June 22, 2014 at 9:18 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1118) said:

IMHO the filtering is optional. It really is just there to try and reduce the
effect of multiple testing. I would check the counts on your differentially
expressed genes and make sure that its not all low count genes, as the calls
on these as DEG might be artifacts and thus are sensitive to these filtering
steps, but I’m not really sure?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1118#respond)

  14. ![](./DESeq2_files/6333281f3f61117325467d71fa1081e0)Stefan on [July 17, 2014 at 6:06 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1345) said:

First at all I want to give you a BIG thank you! I already used your DESeq
tutorial and both of them are very helpful.

When I was using you script, I stumbled accross two errors/problems:

First:  
rownames(mat) <\- colnames(mat) <\- with(colData(dds), paste(condition, type,
sep=" : "))  
produced this error:  
Error in eval(expr, envir, enclos) : object 'type' not found

Second:  
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],  
pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))  
produced this error:  
Error in which(showInPlot) :  
error in evaluating the argument 'x' in selecting a method for function
'which': Error: object 'showInPlot' not found

I used the dataset you provided. Could you help me to figure out, what the
problem is?  
I would appreciate your help!  
Stefan

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1345#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [July 19, 2014 at 2:24 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1370) said:

Thanks for the kind words. Are you using the DESeqDataSetFromHTSeqCount
function to import the table like is shown in the blog? I’m not sure about the
second error (sorry), my guess is your using an newer version or R/DESeq2 and
the function has changed (see my sessionInfo below). I’m a bit rusty on this
(-:

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1370#respond)

      - ![](./DESeq2_files/6333281f3f61117325467d71fa1081e0\(1\))Stefan on [July 21, 2014 at 7:33 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1394) said:

Yes I am using the DESeqDataSetFromHTSeqCount function to import the table. I
downloaded your files and used your script as shown in your blog.

This would be my sessioninfo:  
> sessionInfo()  
R version 3.0.2 (2013-09-25)  
Platform: x86_64-pc-linux-gnu (64-bit)

locale:  
[1] LC_CTYPE=en_US LC_NUMERIC=C LC_TIME=en_US  
[4] LC_COLLATE=en_US LC_MONETARY=en_US LC_MESSAGES=en_US  
[7] LC_PAPER=en_US LC_NAME=C LC_ADDRESS=C  
[10] LC_TELEPHONE=C LC_MEASUREMENT=en_US LC_IDENTIFICATION=C

attached base packages:  
[1] parallel stats graphics grDevices utils datasets methods  
[8] base

other attached packages:  
[1] gplots_2.14.1 RColorBrewer_1.0-5 vsn_3.30.0  
[4] Biobase_2.22.0 DESeq2_1.2.10 RcppArmadillo_0.4.320.0  
[7] Rcpp_0.11.2 GenomicRanges_1.14.4 XVector_0.2.0  
[10] IRanges_1.20.7 BiocGenerics_0.8.0

loaded via a namespace (and not attached):  
[1] affy_1.40.0 affyio_1.30.0 annotate_1.40.1  
[4] AnnotationDbi_1.24.0 BiocInstaller_1.12.1 bitops_1.0-6  
[7] caTools_1.17 DBI_0.2-7 gdata_2.13.3  
[10] genefilter_1.44.0 grid_3.0.2 gtools_3.4.1  
[13] KernSmooth_2.23-12 lattice_0.20-29 limma_3.18.13  
[16] locfit_1.5-9.1 preprocessCore_1.24.0 RSQLite_0.11.4  
[19] splines_3.0.2 stats4_3.0.2 survival_2.37-7  
[22] XML_3.98-1.1 xtable_1.7-3 zlibbioc_1.8.0

  15. ![](./DESeq2_files/6333281f3f61117325467d71fa1081e0)Stefan on [July 23, 2014 at 3:20 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1402) said:

I had these two problems here, right?:

First:  
rownames(mat) <\- colnames(mat) <\- with(colData(dds), paste(condition, type,
sep=" : "))  
produced this error:  
Error in eval(expr, envir, enclos) : object 'type' not found

Second:  
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],  
pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))  
produced this error:  
Error in which(showInPlot) :  
error in evaluating the argument 'x' in selecting a method for function
'which': Error: object 'showInPlot' not found

I am pretty sure you must have gotten the same error, when you ran you script.
For the first error:  
paste(condition, type, sep=" : ") is not working, because at least the type
object is not defined. Otherwise, you figure ([http://dwheelerau.files.wordpre
ss.com/2013/12/deseq2_heatmaps_samplebysample.png?w=584](https://web.archive.o
rg/web/20160328234446/http://dwheelerau.files.wordpress.com/2013/12/deseq2_hea
tmaps_samplebysample.png?w=584)) would have labels seperated by " : " (like
figure 5 in the orignial manual: [http://www.bioconductor.org/packages/2.13/bi
oc/vignettes/DESeq2/inst/doc/DESeq2.pdf](https://web.archive.org/web/201603282
34446/http://www.bioconductor.org/packages/2.13/bioc/vignettes/DESeq2/inst/doc
/DESeq2.pdf))

For the second error:  
You forgot a few lines:  
resFilt <\- res[use & !is.na(res$pvalue),]  
orderInPlot <\- order(resFilt$pvalue)  
showInPlot <\- (resFilt$pvalue[orderInPlot] <= 0.08)  
alpha <\- 0.1

before you can run:  
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],  
pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))

That is the reason, why showInPlot and orderInPlot are undefined…

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1402#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [July 25, 2014 at 8:31 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1411) said:

Hi Stefan, sorry about the errors, thanks for posting some fixes! It really is
a total pain in the arse moving between R and wordpress, adding figures, code,
and rambling on (excuses excuses), so mistakes are bound to pop up. I try to
do these about what I’m doing at work, and I’m not really doing much RNAseq at
the moment (mostly teaching/paperwork). Hopefully I can find some time to
update this next time I get the RNAseq bug (or feel guilty).

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1411#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [July 25, 2014 at 10:25 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1419) said:

I felt guilty and fixed this. Theres one slight issue with the histogram plot
(yellow and blue), not quite sure what I’ve done there. Thanks for the
feedback (and GO K-state wildcats).

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1419#respond)

  16. ![](./DESeq2_files/6b122ea9fe589383445d983e44d1ba8c.jpeg)Anil Kumar on [August 20, 2014 at 2:03 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1647) said:

Hi Sir,

Thank you very much for the descriptive tutorial. I have two doubts. Kindly
help me.

1.In my case, there are 2 control replicates and 2 replicates for each
treatment (cold stress 5 days -> treatment 1 &cold stress 10 days-> treatment
2). So can i take all the different treatment as just treated sample asbelow,

SampleCondition<-c("untreated","untreated","treated","treated","treated","trea
ted")

2\. In the exaple which you shown, sample conditions are specified by
treatment first followed by untreatment. However, in the manual, it is
instructed like this "control" or untreated" level as the first element ("base
level"), so  
that the log2 fold changes produced by default will be the expected comparison
against the base level". Which means we should take control as first
condition?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1647#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [August 20, 2014 at 8:46 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1650) said:

1) I don’t think it would be correct because DESeq will treat your different
treatments as replicates of one treatment rather than separate treatments. You
would need to setup pairwise comparisons, cont v cold stress 5 days and cont V
(cold stress 10 days). Cufflinks has a timeseries options if you wanted to do
what you outlined.

2) Yep, I agree, but see this line
“colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition,
levels=c("untreated","treated"))" I specify untreated first so that the log
values make sense, as you say. I also mentioned this in the text "The levels
in colData are important because they are used in the log calculations; it
makes sense to set untreated or control first so that the direction of the
logs fold changes doesn’t confuse everyone (typically we do comparisons to the
control)! Now for the guts of the DEseq2 analysis.". I think I got it right?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1650#respond)

  17. ![](./DESeq2_files/6b122ea9fe589383445d983e44d1ba8c.jpeg)Anil Kumar on [August 22, 2014 at 2:09 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1665) said:

Hi Sir,  
I wanted to withdraw my second doubt soon after posting it. As you mentioned,
it was clearly indicated in the code as well as in text. Sorry for that..!!

As suggestion by SEQanswers members, Cufflinks and Cuffdiff consider replicate
pair and will not give correct result if am interested in knowing the inter-
replicate variation (actually i need to analyze inter-replicate variation
also) .

As you suggested, i will proceed by splitting my samples as ctrl+cold 5 days
and ctrl+cold 10 days separately. But still have a strong doubt whether such
‘splitting up’ processes is statistically significant when working on
differential expression analysis.

I would like to know your suggestion once again..If that is ok, then i will
convert my count matrix( generated by Featurecount ) to count table (HTSeq)
and proceed as your instructions with spitted data.  
Thanks,  
Anil.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1665#respond)

  18. ![](./DESeq2_files/a457b2c255e5b1e23610ea05f062677b)David on [September 5, 2014 at 3:38 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1833) said:

Hi Dave,

First, like everyone else, I am very grateful for this tutorial. It has helped
me tremendously with analyzing our RNAseq data.

I had a question concerning the heatmap function. Is there away to select for
significantly expressed genes (say those with padj < 0.1) instead of
normalized counts of the 30 highly expressed genes? I was having issues
selecting counts for genes that were differentially expressed. Suggestions?

thank!

library("RColorBrewer")  
library("gplots")  
select <\- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]  
hmcol <\- colorRampPalette(brewer.pal(9, "GnBu"))(100)  
heatmap.2(assay(vsd)[select,], col = hmcol,  
Rowv = FALSE, Colv = FALSE, scale="none",  
dendrogram="none", trace="none", margin=c(10, 6))  
dev.copy(png,"DESeq2_heatmap3")  
dev.off()

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1833#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [September 5, 2014 at 11:21 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1834) said:

Hi There and thanks. I remember i did something like this in DESeq[1]. Not
really near a computer so I can’t play around with it but it might work?? Just
change the select line and keep the rest the same.

select = order(res$padj,decreasing=FALSE)[1:20]  
……

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=1834#respond)

      - ![](./DESeq2_files/a457b2c255e5b1e23610ea05f062677b\(1\))David on [September 8, 2014 at 1:28 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-1857) said:

Hey Dave,

Thanks for the suggestion. I was definitely overthinking the problem. Order by
res$padj works well for me. If any one else is interested, you will need to
add normalization in your heatmap.2 script (scale = ‘row’). e.g.:

select.padj = order(res$padj,decreasing=FALSE)[1:20]  
my_palette <\- colorRampPalette(c("blue",'white','red'))(n=1000)  
par(cex.main=0.8)  
heatmap.2(assay(vsdMF)[select.padj,], col=my_palette,  
scale="row", key=T, keysize=1,symkey=T,density.info="none",  
trace="none",cexCol=0.6, labRow=F)

  19. ![](./DESeq2_files/6c42a4d1bc8ab10499b53bedf375458a)Gonzalo on [January 29, 2015 at 8:12 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3458) said:

thanks for this blog, is quite useful. I am new to DEseq2 and maybe this is
silly question, but I do not understand why meanSdPlot fails…. any ideas to
solve this problem?

library(“vsn”)  
par(mfrow=c(1,3))  
notAllZero 0)  
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))  
meanSdPlot(assay(rld[notAllZero,]))  
meanSdPlot(assay(vsd[notAllZero,]))

Error: could not find function “meanSdPlot”

Thanks in advance

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3458#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [January 30, 2015 at 3:09 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3460) said:

Sounds like the wrong library is being loaded (not sure).  
I get the following when I type meanSdPlot:

>>meanSdPlot  
standardGeneric for “meanSdPlot” defined from package “vsn”

function (x, ranks = TRUE, xlab = ifelse(ranks, “rank(mean)”,  
“mean”), ylab = “sd”, pch = “.”, plot = TRUE, …)  
standardGeneric(“meanSdPlot”)

Methods may be defined for arguments: x, ranks, xlab, ylab, pch, plot  
Use showMethods(“meanSdPlot”) for currently available ones.

My session info from the above included this version of vsn and the
dependencies:  
>>sessionInfo()  
…  
Biobase_2.26.0  
BiocGenerics_0.12.0  
vsn_3.34.0

…

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3460#respond)

      - ![](./DESeq2_files/6c42a4d1bc8ab10499b53bedf375458a\(1\))Gonzalo on [January 30, 2015 at 1:48 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3470) said:

Thanks. I wasnt loading the package “vsn”…

  20. ![](./DESeq2_files/5796b4c580168aadfa66ef0f8fc654cc.jpeg)Anupam Sinha on [April 15, 2015 at 12:43 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3860) said:

Hey Dave,

Thanks for this blog, which is very informative and helpful. I have a query
regarding drawing heatmap. Is it possible to draw a heatmap for about 50 genes
that show greatest change in expression in both directions (i.e upregulation
and downregulation) across samples ? Maybe using by using absolute values of
log2Foldchange in descending order. Thanks for any help…

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3860#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [April 15, 2015 at 8:35 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3862) said:

Hi Sinha, the key line is:  
select <\- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]

As you can see the heatmap is composed of the highest count genes as these
give the most robust results. Changing to greatest fold change you would want
to make sure these don't end up being lots of low count genes with distorted
fold changes. I'm a bit rusty and don't have a dataset to test but play around
with the "select" line until you get what you want, this might be a place to
start:(?)

select <\- order(res$log2FoldChange)[1:30]  
or something like that!

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3862#respond)

  21. ![](./DESeq2_files/508935a56734bed7384605ec8d849e94)Michael Love on [April 22, 2015 at 2:39 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3895) said:

hi David,

Great post. Thanks :)

At some point last year we changed the distance matrix drawing code chunk to:

hc <\- hclust(distsRL)  
heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col =
rev(hmcol), margin=c(13, 13))

This uses the distances directly to order the rows and draw the dendrogram on
top. Otherwise, heatmap.2 will calculate the distances of the distance matrix
to order the rows. (This bug in the vignette was pointed out to me by Jesse
Rowley.)

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3895#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [April 24, 2015 at 11:19 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3905) said:

Hey Michael, thanks for the update and I’m glad I haven’t managed to terribly
misrepresent your great work! And seriously, thanks for DESeq and the great
Documentation that goes with it!

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3905#respond)

  22. ![](./DESeq2_files/a7cc50e9c3699e0ad9ab4f1bd9977ec7)Raheleh Sh on [April 30, 2015 at 7:13 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3923) said:

Dear David,  
Thanks for your post.  
I am trying to make the heatmap only with the most 30 upregulated genes and I
tried this command:  
up 0)  
select.padj <\- order(up$padj,decreasing=FALSE)[1:30]  
hmcol <\- colorRampPalette(brewer.pal(9, "GnBu"))(100)  
pdf("upreg9_heatmap3_vsd_IL22_wt_KO.pdf")  
heatmap.2(assay(vsd)[select.padj,], col = hmcol,  
Rowv = FALSE, Colv = FALSE, scale="column",  
dendrogram="none", trace="none", margin=c(12, 11))  
dev.off()

but, I am getting the heatmap from genes with highest count.  
Would you please help me where is the problem?

Many thanks,  
Rahel

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3923#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [May 6, 2015 at 8:51 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-3936) said:

Sorry Rahel, been really busy. This might help. This is based around res which
is sorted by padj near the top of the code in the post. Using the same logic I
think you could get res to be sorted by +log fold change then use it as
follows.  
#Res is ordered by padj, get top 30  
top30padj<-rownames(res)[1:30]  
s<-which(rownames(dds) %in% top30padj)  
heatmap.2(assay(vsd)[s,], col = hmcol,  
Rowv = FALSE, Colv = FALSE, scale="none",  
dendrogram="none", trace="none", margin=c(10, 8))

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=3936#respond)

  23. ![](./DESeq2_files/c3f2025a4e4180f36772353950a4fdc9)shy on [June 10, 2015 at 8:45 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4101) said:

Hi,I want o know if there are no replicates in my experiments,how can i do ?
Here is my code for replicates. thank you very much!

cds<-DESeqDataSetFromMatrix(data.f3_mf,design=~condition,colData=colData);  
cds<-estimateSizeFactors(cds);  
cds<-estimateDispersions(cds);  
cds<-nbinomWaldTest(cds);  
cds<-DESeq(cds);  
res <\- results(cds)  
resOrdered <\- res[order(res$padj),]  
plotMA(resOrdered);  
addmargins( table( res_sig = resOrdered $padj < .1, res_sig = resOrdered $padj
< .1 ) );

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4101#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [June 10, 2015 at 8:17 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4103) said:

Never done it, but I’d say you just need to change this line:  
sampleCondition<-c("treated","treated","treated","untreated","untreated","untr
eated")

to

sampleCondition<-c("treated","untreated")  
and inport in only the two count files. Remember the value of biological
RNAseq data without replicates is somewhat limited.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4103#respond)

  24. ![](./DESeq2_files/684ba9a9daf0fcc58121bf53159e0bff)[kabasokamfwa](https://web.archive.org/web/20160328234446/http://gravatar.com/kelvinkamfwa) on [June 27, 2015 at 2:52 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4204) said:

Hi, how do I change the padj from the default 0.1 to 0.05

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4204#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [June 28, 2015 at 1:51 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4207) said:

Hi, just pass alpha to the function, most of the functions accept it and the
default is 0.1 ie plotMA(dds,alpha=0.05).

You can see how all the functions work by checking the docs ([http://www.bioco
nductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf](https://web.a
rchive.org/web/20160328234446/http://www.bioconductor.org/packages/release/bio
c/manuals/DESeq2/man/DESeq2.pdf))

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4207#respond)

  25. ![](./DESeq2_files/60151d02d01c6fb2809eedd57fa9debd)jp on [July 19, 2015 at 4:41 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4308) said:

I cant follow your protocol because of this error:  
> ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
directory=directory, design = ~condition)  
Error in DESeqDataSet(se, design = design, ignoreRank) :  
some values in assay are not integers

any help ?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4308#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [July 19, 2015 at 8:56 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4311) said:

It sounds like some of the count values you are using are not integer numbers,
they might be missing values or NaN or are decimals, or perhaps strings (from
a header), I check the count files for anything funny and make sure they
haven’t already been normalised.  
Good luck.  
Cheers

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4311#respond)

  26. ![](./DESeq2_files/2bee8eeaa3f04ac9c28f1b7e6008b521)Little_Jane on [October 27, 2015 at 7:07 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4648) said:

hello, thanks for the good tutorial. I want to know which the method of p
value adjusted is by default. And I also check on the website of DESeq2 and
looked at the command results, which can be used to decide the method we want
to use for adjusting. So, do you know about it?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4648#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [October 27, 2015 at 11:21 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4651) said:

Hi there,  
Its set to BH () for the FDR corrections. If you see the manual ([https://bioc
onductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf](https://web.
archive.org/web/20160328234446/https://bioconductor.org/packages/release/bioc/
manuals/DESeq2/man/DESeq2.pdf)) on pg 35 for the results function there is a
setting pAdjustMethod = “BH”. Based on the doc says you can change it to one
of these:  
p.adjust.methods  
# c(“holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”,  
# “fdr”, “none”)

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4651#respond)

  27. ![](./DESeq2_files/f2b9fe461fc22f8564bcefcc3732dce3)Mayank Kaashyap on [February 18, 2016 at 8:44 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4917) said:

Hi Dave,  
I don’t find your previous version on this. Here you have used amp, quot. Can
I have the previous workflow on DESeq2.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4917#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [February 18, 2016 at 7:46 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4919) said:

Bugger, bloody wordpress – fixed.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4919#respond)

  28. ![](./DESeq2_files/150404ddd38ac7f059593decf7903976)[yifangt](https://web.archive.org/web/20160328234446/http://yifangt.wordpress.com/) on [February 25, 2016 at 7:55 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4935) said:

I am having problem with meanSdPlot() with this error:  
> meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim =
c(0,2.5))  
Error: Unknown parameters: ylim  
> meanSdPlot(assay(rld[notAllZero,]), ylim = c(0,2.5))  
Error: Unknown parameters: ylim  
> meanSdPlot(assay(vsd[notAllZero,]), ylim = c(0,2.5))  
Error: Unknown parameters: ylim  
>

And my sessionInfo() is:  
> sessionInfo()  
R version 3.2.3 (2015-12-10)  
Platform: x86_64-pc-linux-gnu (64-bit)  
Running under: Ubuntu 14.04.4 LTS

locale:  
[1] LC_CTYPE=en_CA.UTF-8 LC_NUMERIC=C  
[3] LC_TIME=en_CA.UTF-8 LC_COLLATE=en_CA.UTF-8  
[5] LC_MONETARY=en_CA.UTF-8 LC_MESSAGES=en_CA.UTF-8  
[7] LC_PAPER=en_CA.UTF-8 LC_NAME=C  
[9] LC_ADDRESS=C LC_TELEPHONE=C  
[11] LC_MEASUREMENT=en_CA.UTF-8 LC_IDENTIFICATION=C

attached base packages:  
[1] parallel stats4 stats graphics grDevices utils datasets  
[8] methods base

other attached packages:  
[1] gplots_2.17.0 RColorBrewer_1.1-2  
[3] vsn_3.38.0 DESeq2_1.10.1  
[5] RcppArmadillo_0.6.500.4.0 Rcpp_0.12.3  
[7] SummarizedExperiment_1.0.2 Biobase_2.30.0  
[9] GenomicRanges_1.22.4 GenomeInfoDb_1.6.3  
[11] IRanges_2.4.7 S4Vectors_0.8.11  
[13] BiocGenerics_0.16.1

loaded via a namespace (and not attached):  
[1] genefilter_1.52.1 gtools_3.5.0 locfit_1.5-9.1  
[4] splines_3.2.3 lattice_0.20-33 colorspace_1.2-6  
[7] survival_2.38-3 XML_3.98-1.3 foreign_0.8-66  
[10] DBI_0.3.1 BiocParallel_1.4.3 affy_1.48.0  
[13] lambda.r_1.1.7 affyio_1.40.0 plyr_1.8.3  
[16] zlibbioc_1.16.0 munsell_0.4.3 gtable_0.1.2  
[19] futile.logger_1.4.1 caTools_1.17.1 labeling_0.3  
[22] latticeExtra_0.6-28 geneplotter_1.48.0 BiocInstaller_1.20.1  
[25] AnnotationDbi_1.32.3 preprocessCore_1.32.0 acepack_1.3-3.3  
[28] KernSmooth_2.23-15 xtable_1.8-2 scales_0.3.0  
[31] limma_3.26.8 gdata_2.17.0 Hmisc_3.17-2  
[34] annotate_1.48.0 XVector_0.10.0 gridExtra_2.0.0  
[37] ggplot2_2.0.0 digest_0.6.9 grid_3.2.3  
[40] bitops_1.0-6 RSQLite_1.0.0 Formula_1.2-1  
[43] cluster_2.0.3 futile.options_1.0.0 rpart_4.1-10  
[46] nnet_7.3-12

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4935#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [February 25, 2016 at 11:02 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4936) said:

That is strange. Your command works on my machine. I checked the manual for
vsn and the meanSdPlot still has that paramater ([http://bioconductor.org/pack
ages/release/bioc/manuals/vsn/man/vsn.pdf](https://web.archive.org/web/2016032
8234446/http://bioconductor.org/packages/release/bioc/manuals/vsn/man/vsn.pdf)
). Perhaps there is a name clash with another library. Try calling it spec
using the package name ie

vsn::meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1), ylim =
c(0,2.5))

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4936#respond)

  29. ![](./DESeq2_files/150404ddd38ac7f059593decf7903976)[yifangt](https://web.archive.org/web/20160328234446/http://yifangt.wordpress.com/) on [February 26, 2016 at 3:15 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4939) said:

Thanks Dave!  
No, calling the function with package name still gave the same error!  
However, it worked if the ylim parameter is removed. Do not know why?!

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4939#respond)

  30. Pingback: [RNAseq pipeline – alignment to DE analysis | Gene](https://web.archive.org/web/20160328234446/https://genehub.wordpress.com/2016/02/29/rnaseq-pipeline-alignment-to-de-analysis/)

  31. ![](./DESeq2_files/493eda420ced32d95a5362678c212e10)[FENG Lei](https://web.archive.org/web/20160328234446/http://genehub.wordpress.com/) on [March 1, 2016 at 9:18 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4955) said:

Reblogged this on [Gene](https://web.archive.org/web/20160328234446/https://ge
nehub.wordpress.com/2016/03/01/how-to-use-deseq2-to-analyse-rnaseq-data-2/)
and commented:  
Good example.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4955#respond)

  32. ![](./DESeq2_files/45598b617fa62796b290d42402eded10)Mayank Kaashyap on [March 7, 2016 at 8:12 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4970) said:

Hi Dave,  
Thanks for such a nice tutorial. May I please know what values are used to
build up a heatmap. Do SIM means correlate to RPKM vlaues. Can I choose few
genes of interest from the table to plot a vsd transformed heatmap?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4970#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 7, 2016 at 7:20 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4972) said:

Hi, I’m not sure I understand your questions. But the heatmap is made from
variance stabilised read counts. I guess the size factor normalised read
counts will correlate with RPKM values, but they are different normalisation
strategies and should not be used interchangeably and RPKM normalised values
should **never** be used with DESeq2. If you want to do a heatmap with just a
few genes, do the variance stabilisation as before ie
vsd<-varianceStabilizingTransformation(dds, blind=TRUE), then you need to
collect the 'index' (probably not R terminology) of the genes you are
interested in and assign them to the value "select" in the tutorial as an
array. If you want to see all the index for the available genes use
rownames(res), use the numbers for each of your genes in your "select"
variable.

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4972#respond)

      - ![](./DESeq2_files/45598b617fa62796b290d42402eded10\(1\))Mayank Kaashyap on [March 8, 2016 at 1:37 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4973) said:

Thanks, Dave!!

  33. ![](./DESeq2_files/45598b617fa62796b290d42402eded10)Mayank on [March 15, 2016 at 1:03 am](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4984) said:

Hi Dave,  
Is there a need to calculate row z-scores?

[Reply
↓](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17
/how-to-use-deseq2-to-analyse-rnaseq-data/?replytocom=4984#respond)

    - ![](./DESeq2_files/1d005972b645b19ff25c553f21888bf5)[Dave Wheeler](https://web.archive.org/web/20160328234446/http://dwheelerau.wordpress.com/) on [March 15, 2016 at 7:00 pm](https://web.archive.org/web/20160328234446/http://dwheelerau.com/2014/02/17/how-to-use-deseq2-to-analyse-rnaseq-data/#comment-4986) said:

I think you are thinking microarrays? The is no need to calculate z-scores.
DESeq normalises the counts and calculates ajusted p-values for the
differential expression calls.  
Cheers

