import PySimpleGUI as sg

sg.theme('Dark Blue 3')

def layout1():
    return [[sg.Text('This is the 1st WINDOW')],
            [sg.Button('Exit1')]]

def layout2():
    return [[sg.Text('This is the 2nd WINDOW')],
            [sg.Button('Exit2')]]

def layoutInit():
    return [[sg.Text('Main WINDOW'),
           sg.Button('Exit')],
          [sg.Button('1st Window'),
           sg.Button('2nd Window')]]


def function_1(values):
    print("caca")

##################Call##################
call_1 = {'1st Window':layout1()}



#########################Main######################
font = ("Helvetica", 16)
window0 = sg.Window('Window Title', layoutInit(), location=(800,200))


active_1, active_2, active_3 = True, False, False
event_1, event_2, event_3 = None, None, None


while True:
    if active_1:
        event_1, values_1 = window0.read(timeout=50)
    if active_2:
        event_2, values_2 = window_2.read(timeout=50)
    if active_3:
        event_3, values_3 = window_3.read(timeout=50)
    # statements for break
    if event_1 in call_1:
        call_1[event_1](values_1)
    # if event_2 in call_2:
    #     call_2[event_2](values1)
    # if event_3 in call_3:
    #     call_3[event_3](values3)
# statements for window close
