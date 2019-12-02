with open(snakemake.input[0]) as reportfile:
    with open(snakemake.output[0], "w") as newreportfile:
        newReport = ""
        for line in reportfile:
            if line.startswith("body {"):
                newReport = "p.cmmd{ text-align: left; padding: 10px; border-style:solid; border-color:#99AAC7; max-width:95%; margin-left:2%; height: auto; background-color: #DEDEFF; word-wrap:normal; }\n commd{text-align: left;} \n span.red{color:red;}\nspan.green{color:#008800;}\n"
                newReport += line
            elif line.startswith("div#metadata {"):
                newReport = "div.document p.cmmd{ text-align: left;}\n"
                newReport += line
            elif "class=\"commd\"" in line:
                newReport = line.replace("<p>","<p class=\"cmmd\">",1)
            else :
                newReport = line
            newreportfile.write(newReport)
        reportfile.close()
        newreportfile.close()
#"p.cmmd{ align: left; padding: 10px; border-style:solid; border-color:#287EC7; max-width:95%; margin-left:2%; height: auto; background-color: #DDDDFF; word-wrap:normal; text-align: left; }"
#<p class="cmmd"><span class="commd">
