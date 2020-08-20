import os
import shutil


def createFolderSave() :

    path = os.getcwd()
    path_final = path+"/savepng_movie"
    os.mkdir(path_final)

    return path_final

def deleteFolder(folder):
    shutil.rmtree(folder)

def savePNG(fig,state,temporary) :

    etat=str(state)
    if temporary == 1 :
        fig.savefig("savepng_movie/etat_"+etat+".png")
    if temporary == 2 :
        fig.savefig("etat"+etat+".png")
