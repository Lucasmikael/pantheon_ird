import PySimpleGUI as sg

sg.theme('Dark Blue 3')

def layout1():
    return [[sg.Text('This is the 1st WINDOW', font=font)],
            [sg.Button('Exit1', font=font)]]

def layout2():
    return [[sg.Text('This is the 2nd WINDOW', font=font)],
            [sg.Button('Exit2', font=font)]]

font = ("Helvetica", 16)
layout = [[sg.Text('Main WINDOW', font=font),
           sg.Button('Exit', font=font)],
          [sg.Button('1st Window', font=font),
           sg.Button('2nd Window', font=font)]]

window0 = sg.Window('Window Title', layout, location=(800,200))
window  = [window0, None, None]
active  = [True, False, False]
event   = [None, None, None]
values  = [None, None, None]

while True:
    for i in range(3):
        if active[i]:
            event[i], values[i] = window[i].read(timeout=50)
            if event[i] != sg.TIMEOUT_KEY:
                print(f'Window {i} event:{event[i]}, values:{values[i]}')
            if event[i] in (sg.WIN_CLOSED, 'Exit', 'Exit1', 'Exit2'):
                if i == 0:
                    for j in range(2,-1, -1):
                        if active[j]:
                            active[j] = False
                            window[j].close()
                    window = None
                    exit()
                else:
                    active[i] = False
                    window[i].close()
            elif event[i] == '1st Window':
                if not active[1]:
                    active[1] = True
                    window[1] = sg.Window("1st Window", layout1(),
                        location=(800, 400), finalize=True)
            elif event[i] == '2nd Window':
                if not active[2]:
                    active[2] = True
                    window[2] = sg.Window("2nd Window", layout2(),
                        location=(800, 600), finalize=True)
