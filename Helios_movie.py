from moviepy.editor import *

def createVideo(folder, file_movie_name):

    list_img = []
    for img in os.listdir(folder):
        list_img.append(folder+"/"+img)
        list_img.sort()
    clips = [ImageClip(m).set_duration(3)
             for m in list_img]

    concat_clip = concatenate_videoclips(clips, method="compose")
    concat_clip.write_videofile(file_movie_name, fps=24)
