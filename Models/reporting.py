import subprocess
import webbrowser

def generatereport(nameproject):
    
    cmd = ["/usr/bin/Rscript" , "Models/generatereports.R " + nameproject]
    subprocess.run(cmd, shell=False, check=True)


def openreportbrowser(nameproject):
    
    pathreport = "Results/Reports/Report_overlapping_" + nameproject + ".html"
    try:
        cmd = ["open", pathreport]
        subprocess.run(cmd, shell=False, check=True)   
    except:
        try:
            webbrowser.open_new(pathreport)
        except:
            print("It was not possible to display automatically the html \
                report.")
