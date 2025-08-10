import numpy as np
import math as mt  # importing libraries
from scipy.linalg import eigh
import cmath as cm
from sympy import symbols, Eq, solve
from math import *
from sympy import *
import sympy as sy
from tkinter.ttk import Button
from tkinter import *
import tkinter as tk
from tkinter import ttk

w = Tk()
tabcontrol = ttk.Notebook(w)

directtab = ttk.Frame(tabcontrol)
rayleightab = ttk.Frame(tabcontrol)
duhameltab = ttk.Frame(tabcontrol)

tabcontrol.add(directtab, text="direct method")
tabcontrol.add(rayleightab, text="rayleigh damping")
tabcontrol.add(duhameltab, text="Duhamel")
tabcontrol.pack(expand=1, fill="both")


def frame1(a, b):
    frame = Frame(directtab, width="480", height="590", background="gray35", bd=2, relief="raised").place(x=a, y=b)


frame = frame1(10, 55)  # for gui(design)

frame3 = Frame(duhameltab, width="900", height="590", background="gray35", bd=2, relief="raised").place(x=55, y=70)


# calculating needed coeficients that are used to find structural response
# A, B, C, D
def coe_C(p0, k, beta, ksi, force_type):
    bet1 = (1 - (beta ** 2))
    bet2 = (2 * ksi * beta)
    pk = (p0 / k)
    if force_type == "sin":
        C = -pk * ((bet2) / ((bet1 ** 2) + (bet2 ** 2)))
    elif force_type == "cosin":
        C = pk * ((bet1) / ((bet1 ** 2) + (bet2 ** 2)))
    else:
        print("  ")
    return C


def coe_D(p0, k, beta, ksi, force_type):
    bet1 = (1 - (beta ** 2))
    bet2 = (2 * ksi * beta)
    pk = (p0 / k)
    if force_type == "sin":
        D = pk * ((bet1) / ((bet1 ** 2) + (bet2 ** 2)))
    elif force_type == "cosin":
        D = -pk * ((bet2) / ((bet1 ** 2) + (bet2 ** 2)))
    else:
        print("  ")
    return D


def coe_A(C, x0):
    A = -C + x0
    return A


def coe_B(xdot0, ksi, wn, A, w, D, wbar):
    B = (xdot0 + (ksi * wn * A) - (w * D)) / wbar
    return B


# forced vibration response
def xt(ksi, wn, t, A, wbar, B, C, D, w):
    e = mt.e
    dt = (e ** (-ksi * wn * t)) * ((A * cm.cos(wbar * t)) + (B * cm.sin(wbar * t))) + (
            (C * cm.cos(w * t)) + (D * cm.sin(w * t)))
    return dt


# now functions for free vibration ...
# B coe for free vibration after turning the engine off

def B_free(xd0, ksi, wn, x0, wbar):
    B = (xd0 + (ksi * wn * x0)) / wbar
    return B


def XT_free(ksi, wn, t, A, wbar, B):
    e = mt.e
    dt = (e ** (-ksi * wn * t)) * ((A * cm.cos(wbar * t)) + (B * cm.sin(wbar * t)))
    return dt


# derivative of forced vibration relation
def XdotT(ksi, wn, T, A, wbar, B, D, C, w):
    e = mt.e
    t = symbols("t")
    dt = (e ** (-ksi * wn * t)) * ((A * sy.cos(wbar * t)) + (B * sy.sin(wbar * t))) + (
            (C * sy.cos(w * t)) + (D * sy.sin(w * t)))
    dt_diff = lambdify(t, dt)
    xdot = dt_diff(T)
    return xdot


def Atriangle(p, w, dt, tfinal, ksi, wn, wbar):
    force_type = forcetype.get()
    f = 0
    if force_type == "sin":
        steps = tfinal / dt
        steps = int(steps)
        for i in range(steps + 1):
            f += (p * sin(w * ((i) * dt))) * (mt.e ** (ksi * wn * dt * (i))) * cos(wbar * dt * (i))
    elif force_type == "cosin":
        steps = tfinal / dt
        steps = int(steps)
        for i in range(steps + 1):
            f += (p * cos(w * ((i) * dt))) * (mt.e ** (ksi * wn * dt * (i))) * cos(wbar * dt * (i))

    AT = f * dt
    return AT


def Btriangle(p, w, dt, tfinal, ksi, wn, wbar):
    force_type = forcetype.get()
    f = 0
    if force_type == "sin":
        steps = tfinal / dt
        steps = int(steps)
        for i in range(steps + 1):
            f += p * mt.sin(w * ((i) * dt)) * (mt.e ** (ksi * wn * dt * (i))) * sin(wbar * dt * (i))
    elif force_type == "cosin":
        steps = tfinal / dt
        steps = int(steps)
        for i in range(steps + 1):
            f += p * mt.cos(w * ((i) * dt)) * (mt.e ** (ksi * wn * dt * (i))) * sin(wbar * dt * (i))

    BT = f * dt
    return BT


def duhammel(ksi, wn, t, m, wbar, AT, BT):
    wn = float(wn)
    ksi = float(ksi)
    t = float(t)
    XT = ((mt.e ** (-ksi * wn * t)) / (m * wbar))*(AT * sin(wbar * t) - BT * cos(wbar * t))
    return XT


m1 = tk.StringVar()
m2 = tk.StringVar()
k1 = tk.StringVar()
k2 = tk.StringVar()
ksi1 = tk.StringVar()
ksi2 = tk.StringVar()
omega1 = tk.StringVar()
omega2 = tk.StringVar()
p01 = tk.StringVar()
p02 = tk.StringVar()
turnof = tk.StringVar()
cutt = tk.StringVar()
x01 = tk.StringVar()
xd01 = tk.StringVar()
x02 = tk.StringVar()
xd02 = tk.StringVar()
forcetype = tk.StringVar()
forcetype.set("cosin")
delta_t = tk.StringVar()

# here inputs or datas will be taken from user

M1 = Entry(directtab, textvariable=m1).place(x="80", y="65")
M2 = Entry(directtab, textvariable=m2).place(x="80", y="90")
K1 = Entry(directtab, textvariable=k1).place(x="80", y="140")
K2 = Entry(directtab, textvariable=k2).place(x="80", y="165")
KSI1 = Entry(directtab, textvariable=ksi1).place(x="80", y="215")
KSI2 = Entry(directtab, textvariable=ksi2).place(x="80", y="240")
OMEGA1 = Entry(directtab, textvariable=omega1).place(x="80", y="290")
OMEGA2 = Entry(directtab, textvariable=omega2).place(x="80", y="315")
P01 = Entry(directtab, textvariable=p01).place(x="270", y="65")
P02 = Entry(directtab, textvariable=p02).place(x="270", y="90")
TURNOF = Entry(directtab, textvariable=turnof).place(x="80", y="365")
CUTT = Entry(directtab, textvariable=cutt, state="disabled", background="white", disabledbackground="gray70")
X01 = Entry(directtab, textvariable=x01).place(x="270", y="215")
XD01 = Entry(directtab, textvariable=xd01).place(x="270", y="240")
X02 = Entry(directtab, textvariable=x02).place(x="270", y="290")
XD02 = Entry(directtab, textvariable=xd02).place(x="270", y="315")
force_type = OptionMenu(directtab, forcetype, "cosin", "sin").place(x=300, y="135")
CUTT.place(x="80", y="390")
X02 = Entry(duhameltab, textvariable=delta_t).place(x="120", y="110")


def activate_entry():
    CUTT.configure(state="normal")


# '''the main function is the code below starts after (def) command and its for structural __
# response by dirext and exact method'''

def dynamic_response():
    M1 = float(m1.get())
    M2 = float(m2.get())
    K1 = float(k1.get())
    K2 = float(k2.get())
    KSI1 = float(ksi1.get())
    KSI2 = float(ksi2.get())
    OMEGA1 = float(omega1.get())
    OMEGA2 = float(omega2.get())
    P01 = float(p01.get())
    P02 = float(p02.get())
    TURNOF = float(turnof.get())
    CUTT = cutt.get()
    X01 = float(x01.get())
    XD01 = float(xd01.get())
    X02 = float(x02.get())
    XD02 = float(xd02.get())

    force_type = forcetype.get()

    # boundary conditions
    BCx0 = [X01, X02]  # distance BC
    BCxd0 = [XD01, XD02]  # speed    BC
    t = TURNOF
    Pt = [P01 * mt.cos(OMEGA1 * t), P02 * mt.cos(OMEGA2 * t)]

    M = np.array([[M1, 0], [0, M2]])
    K = np.array([[K1 + K2, -K2], [-K2, K2]])
    # calculating natural frequencies and mode shape
    V, D = eigh(K, M)
    Wnn = V
    Wn = np.sqrt(V)
    Wn1 = mt.sqrt(Wnn[0])
    Wn2 = mt.sqrt(Wnn[1])

    # making the first row ("1") this is node shape

    D[1] = D[1] * (1 / D[0])
    D[0] = np.ones_like(D[0])

    fi = D
    fit = fi.T

    dot1 = np.dot(fit, M)
    dot2 = np.dot(dot1, fi)
    Mhat = dot2
    Mhat[0][1], Mhat[1][0] = 0, 0

    dot3 = np.dot(fit, K)
    dot4 = np.dot(dot3, fi)
    Khat = dot4
    Khat[0][1], Khat[1][0] = 0, 0
    print("mhat", Mhat)
    print("khat", Khat)
    Chat = np.array([[2 * KSI1 * Wn1 * Mhat[0][0], 0], [0, 2 * KSI2 * Wn2 * Mhat[1][1]]])

    Phat = np.dot(fit, Pt)

    # decoupled equation solve ....
    wn1 = mt.sqrt((Khat[0][0]) / (Mhat[0][0]))
    print("wn1", wn1)
    wn2 = mt.sqrt((Khat[1][1]) / (Mhat[1][1]))
    print("wn2", wn2)
    wbar1 = wn1 * (cm.sqrt(1 - (KSI1 ** 2)))
    print("wbar1", wbar1)
    wbar2 = wn2 * (cm.sqrt(1 - (KSI2 ** 2)))
    print("wbar2", wbar2)
    fiinv = np.linalg.inv(fi)
    BCq0 = np.dot(BCx0, fiinv)
    BCqd0 = np.dot(BCxd0, fiinv)
    phat011 = P01
    phat012 = (fit[0][1]) * P02
    phat021 = P01
    phat022 = (fit[1][1]) * P02
    # the general format of the forces in decoupled equations ... =                             <<<<<<<<<<<
    # p01cos (w1*t)+ phat012 cos(w2*t)   and the next line will be like this....               <<<<<<<<<<<
    # for the first part >>>> forced vibration ...
    C11 = coe_C(p0=phat011, k=Khat[0][0], beta=(OMEGA1 / wn1), ksi=KSI1, force_type=force_type)
    print("C11", C11)
    A11 = coe_A(C=C11, x0=BCq0[0])
    print("A11", A11)
    D11 = coe_D(p0=phat011, k=Khat[0][0], beta=(OMEGA1 / wn1), ksi=KSI1, force_type=force_type)
    print("D11", D11)
    B11 = coe_B(xdot0=BCqd0[0], ksi=KSI1, wn=wn1, A=A11, w=OMEGA1, D=D11, wbar=wbar1)
    print("B11", B11)
    e = mt.e
    q11 = xt(ksi=KSI1, wn=wn1, t=TURNOF, A=A11, wbar=wbar1, B=B11, C=C11, D=D11, w=OMEGA1)
    print("q11", q11)
    # this was q1 acording to the first force in decoupled equations

    C12 = coe_C(p0=phat012, k=Khat[0][0], beta=(OMEGA2 / wn1), ksi=KSI1, force_type=force_type)
    A12 = coe_A(C=C12, x0=BCq0[0])
    D12 = coe_D(p0=phat012, k=Khat[0][0], beta=(OMEGA2 / wn1), ksi=KSI1, force_type=force_type)
    B12 = coe_B(xdot0=BCqd0[0], ksi=KSI1, wn=wn1, A=A12, w=OMEGA2, D=D12, wbar=wbar1)
    e = mt.e
    q12 = xt(ksi=KSI1, wn=wn1, t=TURNOF, A=A12, wbar=wbar1, B=B12, C=C12, D=D12, w=OMEGA2)
    print("q12", q12)
    # and this was for first decoupled equation according to second phat force

    q1 = (q11 + q12)
    print("q1", q1)
    C21 = coe_C(p0=phat011, k=Khat[1][1], beta=(OMEGA1 / wn2), ksi=KSI2, force_type=force_type)
    print("C21", C21)
    A21 = coe_A(C=C21, x0=BCq0[1])
    print("A21", A21)
    D21 = coe_D(p0=phat011, k=Khat[1][1], beta=(OMEGA1 / wn2), ksi=KSI2, force_type=force_type)
    print("D21", D21)
    B21 = coe_B(xdot0=BCqd0[1], ksi=KSI2, wn=wn2, A=A21, w=OMEGA1, D=D21, wbar=wbar2)
    print("B21", B21)
    e = mt.e
    q21 = xt(ksi=KSI2, wn=wn2, t=TURNOF, A=A21, wbar=wbar2, B=B21, C=C21, D=D21, w=OMEGA1)
    print("q21", q21)
    # this was q2 acording to the first force in decoupled equations

    C22 = coe_C(p0=phat022, k=Khat[1][1], beta=(OMEGA2 / wn2), ksi=KSI2, force_type=force_type)
    A22 = coe_A(C=C22, x0=BCq0[1])
    D22 = coe_D(p0=phat022, k=Khat[1][1], beta=(OMEGA2 / wn2), ksi=KSI2, force_type=force_type)
    B22 = coe_B(xdot0=BCqd0[1], ksi=KSI2, wn=wn2, A=A22, w=OMEGA2, D=D22, wbar=wbar2)
    e = mt.e
    q22 = xt(ksi=KSI2, wn=wn2, t=TURNOF, A=A22, wbar=wbar2, B=B22, C=C22, D=D22, w=OMEGA2)

    q2 = q21 + q22
    print("q2", q2)
    q = [q1, q2]
    print("q", q)
    XT = np.dot(fi, q)  # this is XT in forced vibration time....
    print("XT", XT)
    if CUTT != '':
        CUTT = float(CUTT)
        XdotT11 = XdotT(ksi=KSI1, wn=wn1, T=TURNOF, A=A11, wbar=wbar1, B=B11, C=C11, D=D11, w=OMEGA1)
        XdotT12 = XdotT(ksi=KSI1, wn=wn1, T=TURNOF, A=A12, wbar=wbar1, B=B12, C=C12, D=D12, w=OMEGA2)
        XdotT1 = XdotT11 + XdotT12
        print("xdot1",XdotT1 )
        XdotT21 = XdotT(ksi=KSI2, wn=wn2, T=TURNOF, A=A21, wbar=wbar2, B=B21, C=C21, D=D21, w=OMEGA1)
        XdotT22 = XdotT(ksi=KSI2, wn=wn2, T=TURNOF, A=A22, wbar=wbar2, B=B22, C=C22, D=D22, w=OMEGA2)
        XdotT2 = XdotT21 + XdotT22
        print("xdot2", XdotT2)
        BF11 = B_free(xd0=XdotT1, ksi=KSI1, wn=wn1, x0=q1, wbar=wbar1)
        BF22 = B_free(xd0=XdotT2, ksi=KSI2, wn=wn2, x0=q2, wbar=wbar2)
        print("BF1", BF11)
        print("BF2", BF22)
        AF11 = q1
        AF22 = q2
        print("a1", AF11)
        print("a2", AF22)
        qf1 = XT_free(ksi=KSI1, wn=wn1, t=(CUTT - TURNOF), A=AF11, wbar=wbar1, B=BF11)
        qf2 = XT_free(ksi=KSI2, wn=wn2, t=(CUTT - TURNOF), A=AF22, wbar=wbar2, B=BF22)
        print("qf1", qf1)
        print("qf2", qf2)
        QF = [qf1, qf2]
        print("qf", QF)

        XF = np.dot(fi, QF)
        print("XF", XF)
    else:
        XF = ['', '']

    return XT, XF


def solveduhamel():
    M1 = float(m1.get())
    M2 = float(m2.get())
    K1 = float(k1.get())
    K2 = float(k2.get())
    KSI1 = float(ksi1.get())
    KSI2 = float(ksi2.get())
    OMEGA1 = float(omega1.get())
    OMEGA2 = float(omega2.get())
    P01 = float(p01.get())
    P02 = float(p02.get())
    TURNOF = float(turnof.get())
    CUTT = cutt.get()
    X01 = float(x01.get())
    XD01 = float(xd01.get())
    X02 = float(x02.get())
    XD02 = float(xd02.get())
    Delta_t = float(delta_t.get())

    force_type = forcetype.get()

    # boundary conditions
    BCx0 = [X01, X02]  # distance BC
    BCxd0 = [XD01, XD02]  # speed    BC
    t = TURNOF
    Pt = [P01 * mt.cos(OMEGA1 * t), P02 * mt.cos(OMEGA2 * t)]

    M = np.array([[M1, 0], [0, M2]])
    K = np.array([[K1 + K2, -K2], [-K2, K2]])
    # calculating natural frequencies and mode shape
    V, D = eigh(K, M)
    Wnn = V
    Wn = np.sqrt(V)
    Wn1 = mt.sqrt(Wnn[0])
    Wn2 = mt.sqrt(Wnn[1])

    # making the first row ("1") this is node shape

    D[1] = D[1] * (1 / D[0])
    D[0] = np.ones_like(D[0])

    fi = D
    fit = fi.T

    dot1 = np.dot(fit, M)
    dot2 = np.dot(dot1, fi)
    Mhat = dot2
    Mhat[0][1], Mhat[1][0] = 0, 0

    dot3 = np.dot(fit, K)
    dot4 = np.dot(dot3, fi)
    Khat = dot4
    Khat[0][1], Khat[1][0] = 0, 0

    Chat = np.array([[2 * KSI1 * Wn1 * Mhat[0][0], 0], [0, 2 * KSI2 * Wn2 * Mhat[1][1]]])

    Phat = np.dot(fit, Pt)

    # decoupled equation solve ....
    wn1 = mt.sqrt((Khat[0][0]) / (Mhat[0][0]))
    wn2 = mt.sqrt((Khat[1][1]) / (Mhat[1][1]))
    wbar1 = wn1 * (cm.sqrt(1 - (KSI1 ** 2)))
    wbar2 = wn2 * (cm.sqrt(1 - (KSI2 ** 2)))
    fiinv = np.linalg.inv(fi)
    BCq0 = np.dot(BCx0, fiinv)
    BCqd0 = np.dot(BCxd0, fiinv)
    phat011 = P01
    phat012 = (fit[0][1]) * P02
    phat021 = P01
    phat022 = (fit[1][1]) * P02
    A11 = Atriangle(p=phat011, w=OMEGA1, dt=Delta_t, tfinal=TURNOF, ksi=KSI1, wn=wn1, wbar=wbar1)
    B11 = Btriangle(p=phat011, w=OMEGA1, dt=Delta_t, tfinal=TURNOF, ksi=KSI1, wn=wn1, wbar=wbar1)
    q11 = duhammel(ksi=KSI1, wn=wn1, t=TURNOF, m=Mhat[0][0], wbar=wbar1, AT=A11, BT=B11)

    A12 = Atriangle(p=phat012, w=OMEGA2, dt=Delta_t, tfinal=TURNOF, ksi=KSI1, wn=wn1, wbar=wbar1)
    B12 = Btriangle(p=phat012, w=OMEGA2, dt=Delta_t, tfinal=TURNOF, ksi=KSI1, wn=wn1, wbar=wbar1)
    q12 = duhammel(ksi=KSI1, wn=wn1, t=TURNOF, m=Mhat[0][0], wbar=wbar1, AT=A12, BT=B12)

    q1 = q11 + q12

    A21 = Atriangle(p=phat011, w=OMEGA1, dt=Delta_t, tfinal=TURNOF, ksi=KSI2, wn=wn2, wbar=wbar2)
    B21 = Btriangle(p=phat011, w=OMEGA1, dt=Delta_t, tfinal=TURNOF, ksi=KSI2, wn=wn2, wbar=wbar2)
    q21 = duhammel(ksi=KSI2, wn=wn2, t=TURNOF, m=Mhat[1][1], wbar=wbar2, AT=A21, BT=B21)

    A22 = Atriangle(p=phat022, w=OMEGA2, dt=Delta_t, tfinal=TURNOF, ksi=KSI2, wn=wn2, wbar=wbar2)
    B22 = Btriangle(p=phat022, w=OMEGA2, dt=Delta_t, tfinal=TURNOF, ksi=KSI2, wn=wn2, wbar=wbar2)
    q22 = duhammel(ksi=KSI2, wn=wn2, t=TURNOF, m=Mhat[1][1], wbar=wbar2, AT=A22, BT=B22)
    q2 = q21 + q22
    q = [q1, q2]
    XT = np.dot(fi, q)

    # if CUTT != '':
    #     qf11 = duhammel(ksi=KSI1, wn=wn1, t=CUTT, m=Mhat[0][0], wbar=wbar1, AT=A11, BT=B11)
    #     qf12 = duhammel(ksi=KSI1, wn=wn1, t=CUTT, m=Mhat[0][0], wbar=wbar1, AT=A12, BT=B12)
    #     qf1 = qf11 + qf12
    #
    #     qf21 = duhammel(ksi=KSI2, wn=wn2, t=CUTT, m=Mhat[1][1], wbar=wbar2, AT=A21, BT=B21)
    #     qf22 = duhammel(ksi=KSI2, wn=wn2, t=CUTT, m=Mhat[1][1], wbar=wbar2, AT=A22, BT=B22)
    #     qf2 = qf21 + qf22
    #     qf = [qf1, qf2]
    #     XTF = np.dot(fi, qf)
    #
    # else:
    #     XTF = ["", ""]

    return XT


def outputDu():
    x = solveduhamel()
    TURNOF = float(turnof.get())
    ex1 = x[0]
    ex2 = x[1]
    story1 = Label(duhameltab, text=f"first story displacement : \n  {x} in time :{TURNOF} s", bg="gray25",
                   fg="ghost white",
                   font="arvo 10 bold italic", relief="sunken").place(x="200", y="200")
    story2 = Label(duhameltab, text=f"second story displacement : \n  {x} in time :{TURNOF} s", bg="gray25",
                   fg="ghost white",
                   font="arvo 10 bold italic", relief="sunken").place(x="200", y="300")


# the function is for rayleigh damping that can be used to get the information such as alpha, beta, C and modal C



def rayleigh_damping():
    M1 = float(m1.get())
    M2 = float(m2.get())
    K1 = float(k1.get())
    K2 = float(k2.get())
    KSI1 = float(ksi1.get())
    KSI2 = float(ksi2.get())
    M = np.array([[M1, 0], [0, M2]])
    K = np.array([[K1 + K2, -K2], [-K2, K2]])
    # calculating natural frequencies and mode shape
    V, D = eigh(K, M)
    Wnn = V
    Wn = np.sqrt(V)
    Wn1 = mt.sqrt(Wnn[0])
    Wn2 = mt.sqrt(Wnn[1])

    D[1] = D[1] * (1 / D[0])
    D[0] = np.ones_like(D[0])

    fi = D
    fit = fi.T

    alpha, beta = symbols("alpha, beta")
    eq1 = Eq((alpha + (beta * (Wn1 ** 2))), (2 * KSI1 * Wn1))
    eq2 = Eq((alpha + (beta * (Wn2 ** 2))), (2 * KSI2 * Wn2))
    alphabeta = solve((eq1, eq2), (alpha, beta))
    alpha = float(alphabeta[alpha])
    beta = float(alphabeta[beta])
    CR = np.dot(alpha, M) + np.dot(beta, K)
    dot1 = np.dot(fit, CR)
    dot2 = np.dot(dot1, fi)
    CHR = dot2
    CHR[0][1] = 0
    CHR[1][0] = 0
    l = Label(rayleightab,
              text=f"alpha ={alpha}\n\nbeta ={beta}\n\n\ndamping matrix = \n{CR}\n\n\nmodal damping matrix =\n{CHR}\n\n\nC = (alpha x M)+(beta x K)",
              background="gray35", bd=2, relief="raised", font="arvo 16 bold italic", fg="white", justify="left")
    l.place(x=10, y=150)


def data1():
    M1 = float(m1.get())
    M2 = float(m2.get())
    K1 = float(k1.get())
    K2 = float(k2.get())
    M = np.array([[M1, 0], [0, M2]])
    K = np.array([[K1 + K2, -K2], [-K2, K2]])
    # calculating natural frequencies and mode shape
    V, D = eigh(K, M)
    Wnn = V
    Wn = np.sqrt(V)
    Wn1 = mt.sqrt(Wnn[0])
    Wn2 = mt.sqrt(Wnn[1])
    D[1] = D[1] * (1 / D[0])
    D[0] = np.ones_like(D[0])
    fi = D
    l = Label(directtab, text=f"natural frequency (Wn1) : \n  {Wn1} ", bg="gray25",
              fg="ghost white",
              font="arvo 10 bold italic", relief="sunken").place(x="520", y="300")
    l2 = Label(directtab, text=f"natural frequency (Wn2) : \n  {Wn2} ", bg="gray25",
               fg="ghost white",
               font="arvo 10 bold italic", relief="sunken").place(x="520", y="350")
    l3 = Label(directtab, text=f" modal shape matrix : \n  {fi} ", bg="gray25",
               fg="ghost white",
               font="arvo 10 bold italic", relief="sunken").place(x="520", y="400")


# this function is for outputs
def output1():
    TURNOF = float(turnof.get())
    x = dynamic_response()
    ex1 = x[0][0]
    ex2 = x[0][1]
    exf1 = x[1][0]
    exf2 = x[1][1]
    story1 = Label(directtab, text=f"first story displacement : \n  {ex1} in time :{TURNOF} s", bg="gray25",
                   fg="ghost white",
                   font="arvo 10 bold italic", relief="sunken").place(x="520", y="200")
    story2 = Label(directtab, text=f"second story displacement : \n {ex2} in time :{TURNOF} s", bg="gray25",
                   fg="ghost white",
                   font="arvo 10 bold italic", relief="sunken").place(x="520", y="250")
    if exf1 != '':
        story1f = Label(directtab, text=f"first story displacement :\n  {exf1} in cut time (free vibration)",
                        bg="gray25",
                        fg="ghost white",
                        font="arvo 10 bold italic", relief="sunken").place(x="520", y="70")
        story2f = Label(directtab, text=f"second story displacement : \n    {exf2} in cut time (free vibration) ",
                        bg="gray25",
                        fg="ghost white",
                        font="arvo 10 bold italic", relief="sunken").place(x="520", y="120")
    else:
        aaa = None


# codes below are for window design

w.geometry("1000x720")
w.title("dynamics of structures/ hossein bagheri ")
l = Label(directtab, text="dynamic response of 2dof---force type:p0sin(wt)\ngithub: hosseinbagheri0110 ",
          fg="ghost white",
          bg="gray20"
          , font="arvo 13 bold italic", bd=3, relief="raised", width="99").place(x=2, y=7)
lr = Label(rayleightab, text="dynamic response of 2dof---force type:p0sin(wt)\ngithub: hosseinbagheri0110 ",
           fg="ghost white",
           bg="gray20"
           , font="arvo 13 bold italic", bd=3, relief="raised", width="99").place(x=2, y=7)
lD = Label(duhameltab, text="dynamic response of 2dof---force type:p0sin(wt)\ngithub: hosseinbagheri0110 ",
           fg="ghost white",
           bg="gray20"
           , font="arvo 13 bold italic", bd=3, relief="raised", width="99").place(x=2, y=7)

l2 = Label(directtab, text="email : hosseinbagheri0110@gmail.com", fg="ghost white",
           bg="gray20"
           , font="comic 12 bold italic", bd=3, relief="raised", width="30").place(x=670, y=660)

w.configure(background="gray45")
w.resizable(width=False, height=False)
# button1 = tk.Button(text="submit", width=25, activeforeground="white", activebackground="black", bg="gray25",
#                     fg="ghost white", font="arvo 16 bold italic", command=w.destroy).place(x=75, y=650)

lablem1 = Label(directtab, text="M1    ", background="gray27", fg="white", font="arvo 10 bold italic", bd=2,
                relief="sunken").place(x="25", y="65")
lablem2 = Label(directtab, text="M2    ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                background="gray27").place(x="25", y="90")
lablek1 = Label(directtab, text="k1    ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                background="gray27").place(x="25", y="140")
lablek2 = Label(directtab, text="k2    ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                background="gray27").place(x="25", y="165")
lableksi1 = Label(directtab, text="ksi 1 ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="25", y="215")
lableksi2 = Label(directtab, text="ksi 2", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="25", y="240")
lablew1 = Label(directtab, text="omega1 ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                background="gray27").place(x="19", y="290")
lablew2 = Label(directtab, text="omega2 ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                background="gray27").place(x="19", y="315")
lablep01 = Label(directtab, text="p(0)1 ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                 background="gray27").place(x="230", y="65")
lablep02 = Label(directtab, text="p(0)2 ", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                 background="gray27").place(x="230", y="90")
lableofft = Label(directtab, text="Time 1", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="15", y="365")
lablecutt = Label(directtab, text="Time 2", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="15", y="390")
lablex01 = Label(directtab, text="x01", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                 background="gray27").place(x="230", y="215")
lablexd01 = Label(directtab, text="x.01", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="230", y="240")
labletype = Label(directtab, text="force type", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="230", y="140")
lablex02 = Label(directtab, text="x.02", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                 background="gray27").place(x="230", y="290")
lablexd02 = Label(directtab, text="x.02", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="230", y="315")
buttonof = Button(directtab, text="engine has turned off and\nthere is free vibration", fg="white",
                  font="arvo 8 bold italic", bd=2, relief="raised",
                  background="gray27", cursor="hand2", command=activate_entry).place(x="220", y="368")

labledelta = Label(duhameltab, text="delta t", fg="white", font="arvo 10 bold italic", bd=2, relief="sunken",
                  background="gray27").place(x="70", y="110")


frame2 = frame1(510, 55)
button1 = tk.Button(directtab, text="submit", width=25, activeforeground="white", activebackground="black", bg="gray25",
                    fg="ghost white", font="arvo 16 bold italic", cursor="hand2",
                    command=lambda: [output1(), data1()]).place(x=75, y=650)

button2 = tk.Button(rayleightab, text="calculate Rayleigh Damping according to direct method tab data", width=60,
                    activeforeground="white", activebackground="black", bg="gray25",
                    fg="ghost white", font="arvo 16 bold italic", cursor="hand2", command=rayleigh_damping).place(x=75,
                                                                                                                  y=60)
button3 = tk.Button(duhameltab, text="calculate X(t) with duhamel", width=60,
                    activeforeground="white", activebackground="black", bg="gray25",
                    fg="ghost white", font="arvo 16 bold italic", cursor="hand2", command=outputDu).place(x=75,
                                                                                                                  y=60)
w.mainloop()
