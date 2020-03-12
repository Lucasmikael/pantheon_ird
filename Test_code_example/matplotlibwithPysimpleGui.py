import PySimpleGUI as sg
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")

"""
    Simultaneous PySimpleGUI Window AND a Matplotlib Interactive Window
    A number of people have requested the ability to run a normal PySimpleGUI window that
    launches a MatplotLib window that is interactive with the usual Matplotlib controls.
    It turns out to be a rather simple thing to do.  The secret is to add parameter block=False to plt.show()
"""

def draw_plot():
    plt.plot([0.1, 0.2, 0.5, 0.7])
    plt.show(block=False)

layout = [[sg.Button('Plot'), sg.Cancel(), sg.Button('Popup')]]

window = sg.Window('Have some Matplotlib....', layout)

while True:
    event, values = window.read()
    if event in (None, 'Cancel'):
        break
    elif event == 'Plot':
        draw_plot()
    elif event == 'Popup':
        sg.popup('Yes, your application is still running')
window.close()

# import PySimpleGUI as sg
# import numpy as np
# import matplotlib
# matplotlib.use('Qt5Agg')
# import matplotlib.pyplot as plt
# from scipy import interpolate
#
#
# def kWindow(valores=None):
#     line = "[sg.In('{0}', size=(20,1), key='-temp{id}-', visible={2}, justification='right'), sg.In('{1}', size=(20,1), key='-k{id}-', visible={2}, justification='right'), sg.Button('+', size=(5,1), key='+{id}', visible={3})]"
#     layout = [[sg.T('x', size=(20,1), justification='center'), sg.T('y', size=(20,1), justification='center')]]
#     layout.append([sg.Button('Plot', size=(20,1)), sg.Button('Cancel', size=(20,1))])
#     for i in range(0, 10):
#         try:
#             if valores is None:
#                 layout.insert(-1, eval(line.format('','','True','True',id=i)))
#                 valores = np.array([['','']])
#             elif i+1 == valores.shape[0]:
#                 layout.insert(-1, eval(line.format(valores[i,0],valores[i,1],'True','True',id=i)))
#             else:
#                 layout.insert(-1, eval(line.format(valores[i,0],valores[i,1],'True','False',id=i)))
#         except: layout.insert(-1, eval(line.format('','','False','False',id=i)))
#     window = sg.Window('Some Name', layout)
#     while True:
#         event, values = window.read()
#         if event in (None, 'Exit', 'Cancel'):
#             window.close()
#             break
#         if "Plot" in event:
#             try:
#                 hr, t = [], []
#                 valores = list(values.values())
#                 for i in range(0, len(valores)-1, 2):
#                     if valores[i]!='' and valores[i+1]!='':
#                         hr.append(valores[i])
#                         t.append(valores[i+1])
#                 params = np.array(list(zip(hr, t)), dtype=float)
#                 plotk(params)
#             except ValueError:
#                 sg.PopupOK('Invalid Input!', title='Warning')
#         if "+" in event:
#             n = int(event[1:])+1
#             if n<9:
#                 window[f'+{n-1}'].update(visible=False)
#                 window[f'-temp{n}-'].update(visible=True)
#                 window[f'-k{n}-'].update(visible=True)
#                 window[f'+{n}'].update(visible=True)
#             elif n==9:
#                 window[f'+{n-1}'].update(visible=False)
#                 window[f'-temp{n}-'].update(visible=True)
#                 window[f'-k{n}-'].update(visible=True)
#
# def plotk(params):
#     T, k = params.T
#     s = interpolate.InterpolatedUnivariateSpline(T, k, k=2)
#     x = np.linspace(T.min(), T.max(), 100)
#     y = s(x)
#     plt.scatter(T, k)
#     plt.plot(x, y)
#     plt.show()
#
#
# mockdata = np.array([[25, 7],[150, 6], [400, 5], [700, 5.5], [1200, 6.5]])
# kWindow(mockdata)
