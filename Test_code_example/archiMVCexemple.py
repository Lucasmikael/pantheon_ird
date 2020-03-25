# from tkinter import Tk, Toplevel, Label, Entry, Button
#
#
# class Transaction:
#     def __init__(self):
#         self.value = 0
#         self.callbacks = {}
#
#     def add_callback(self, func):
#         self.callbacks[func] = None
#
#     def _callbacks(self):
#         for func in self.callbacks:
#             func(self.value)
#
#     def set(self, value):
#         self.value = value
#         self._callbacks()
#
#     def get(self):
#         return self.value
#
#
# class Account:  # The Model
#     def __init__(self):
#         self.transaction = Transaction()
#
#     def deposit(self, value):
#         self.transaction.set(self.transaction.get() + value)
#
#     def withdrawal(self, value):
#         self.transaction.set(self.transaction.get() - value)
#
#
# class BankView(Toplevel):  # View 1
#     def __init__(self, master):
#         Toplevel.__init__(self, master)
#         self.protocol('WM_DELETE_WINDOW', self.master.destroy)
#
#         Label(self, text='Account Balance').pack(side='left')
#         self.balance = Entry(self, width=8)
#         self.balance.pack(side='left')
#
#     def set_balance(self, amount):
#         self.balance.delete(0, 'end')
#         self.balance.insert('end', str(amount))
#
#
# class TellerView(Toplevel):  # View 2
#     def __init__(self, master):
#         Toplevel.__init__(self, master)
#
#         Label(self, text='Amount').pack(side='left')
#         self.amount = Entry(self, width=8)
#         self.amount.pack(side='left')
#
#         self.btn_deposit = Button(self, text='Deposit', width=8)
#         self.btn_deposit.pack(side='left')
#
#         self.btn_withdrawal = Button(self, text='Withdrawal', width=8)
#         self.btn_withdrawal.pack(side='left')
#
#
# class Bank(Tk):  # The Controller
#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.withdraw()
#
#         self.account = Account()
#         self.account.transaction.add_callback(self.update_account)
#
#         self.bank_view = BankView(self)
#         self.bank_view.title('The Bank')
#
#         self.teller_view = TellerView(self.bank_view)
#         self.teller_view.title('The Teller')
#
#         self.teller_view.btn_deposit.config(command=self.make_deposit)
#         self.teller_view.btn_withdrawal.config(command=self.make_withdrawal)
#
#         self.update_account(self.account.transaction.get())
#
#     def make_deposit(self):
#         self.account.deposit(int(self.teller_view.amount.get()))
#
#     def make_withdrawal(self):
#         self.account.withdrawal(int(self.teller_view.amount.get()))
#
#     def update_account(self, amount):
#         self.bank_view.set_balance(amount)
#
#
# if __name__ == '__main__':
#     Bank().mainloop()


import numpy as np
import matplotlib.pyplot as plt

class DraggableRectangle:
    def __init__(self, rect):
        self.rect = rect
        self.press = None

    def connect(self):
        'connect to all the events we need'
        self.cidpress = self.rect.figure.canvas.mpl_connect(
            'button_press_event', self.on_press)
        self.cidrelease = self.rect.figure.canvas.mpl_connect(
            'button_release_event', self.on_release)
        self.cidmotion = self.rect.figure.canvas.mpl_connect(
            'motion_notify_event', self.on_motion)

    def on_press(self, event):
        'on button press we will see if the mouse is over us and store some data'
        if event.inaxes != self.rect.axes: return

        contains, attrd = self.rect.contains(event)
        if not contains: return
        print('event contains', self.rect.xy)
        x0, y0 = self.rect.xy
        self.press = x0, y0, event.xdata, event.ydata

    def on_motion(self, event):
        'on motion we will move the rect if the mouse is over us'
        if self.press is None: return
        if event.inaxes != self.rect.axes: return
        x0, y0, xpress, ypress = self.press
        dx = event.xdata - xpress
        dy = event.ydata - ypress
        #print('x0=%f, xpress=%f, event.xdata=%f, dx=%f, x0+dx=%f' %
        #      (x0, xpress, event.xdata, dx, x0+dx))
        self.rect.set_x(x0+dx)
        self.rect.set_y(y0+dy)

        self.rect.figure.canvas.draw()


    def on_release(self, event):
        'on release we reset the press data'
        self.press = None
        self.rect.figure.canvas.draw()

    def disconnect(self):
        'disconnect all the stored connection ids'
        self.rect.figure.canvas.mpl_disconnect(self.cidpress)
        self.rect.figure.canvas.mpl_disconnect(self.cidrelease)
        self.rect.figure.canvas.mpl_disconnect(self.cidmotion)

fig = plt.figure()
ax = fig.add_subplot(111)
rects = ax.bar(range(10), 20*np.random.rand(10))
drs = []
for rect in rects:
    dr = DraggableRectangle(rect)
    dr.connect()
    drs.append(dr)

plt.show()
