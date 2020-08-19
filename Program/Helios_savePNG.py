import os
import shutil


def createFolderSave() :

    path = os.getcwd()
    path_final = path+"/savepng_movie"
    os.mkdir(path_final)

    return path_final

def deleteFolder(folder):
    shutil.rmtree(folder)

def savePNG(fig,state) :

    etat=str(state)
    fig.savefig("savepng_movie/etat_"+etat+".png")
