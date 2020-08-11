from moviepy.editor import *

def createVideo():
    # clip = VideoClip(make_frame = None, duration=4)  # for custom animations (see below)
    # clip = VideoFileClip("my_video_file.mp4")  # or .avi, .webm, .gif ...
    # clip = ImageSequenceClip(['0 State.png', "1 State.png", "2 State.png"], fps=24)
    # clip = ImageClip(".png")  # or .jpeg, .tiff, ...
    # clip = TextClip("Hello !", font="Amiri-Bold", fontsize=70, color="black")
    # clip = ColorClip(size=(460, 380))

    img = ['0 State.png', '1 State.png', '2 State.png']

    clips = [ImageClip(m).set_duration(3)
             for m in img]

    concat_clip = concatenate_videoclips(clips, method="compose")
    concat_clip.write_videofile("test.mp4", fps=24)


