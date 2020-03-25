# import PySimpleGUI as sg
#
#
# layout_start =  [
#     [sg.Radio('My first Radio!', "RADIO1", default=True),
#     sg.Radio('My second radio!', "RADIO1")]
# ]
#
# window_start = sg.Window("Initialize Network", alpha_channel=0.95, layout=layout_start)
# event_start, values_start = window_start.Read()
#

import sys

if sys.version_info[0] >= 3:
    import PySimpleGUI as sg
else:
    import PySimpleGUI27 as sg

buttons_col = []
for i in range(5):
    buttons_col.append([sg.Button('col {}'.format(i))])

buttons_row = []
for i in range(5):
    buttons_row.append(sg.Button('row {}'.format(i)))

layout = [
            [sg.Text('Your typed chars appear here:'), sg.Text('', key='_OUTPUT_')],
            *buttons_col,
            [*buttons_row],
            [sg.Input(do_not_clear=True, key='_IN_')],
            [sg.Button('Show'), sg.Button('Exit')]
         ]

window = sg.Window('Window Title').Layout(layout)

while True:             # Event Loop
    event, values = window.Read()
    print(event, values)
    if event is None or event == 'Exit':
        break
    if event == 'Show':
        # change the "output" element to be the value of "input" element
        window.FindElement('_OUTPUT_').Update(values['_IN_'])

window.Close()
