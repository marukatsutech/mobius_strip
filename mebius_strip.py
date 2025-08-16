""" Möbius strip """
import numpy as np
from matplotlib.figure import Figure
import matplotlib.animation as animation
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
import tkinter as tk
from tkinter import ttk
from scipy.spatial.transform import Rotation
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d import proj3d

""" Global variables """
cross_lines = []
num_cross_lines = 240
pitch_spiral = 0.1
num_lap = 1.
amplitude_v = 0.5
amplitude_h = 0.5
num_twist = 0.5
wave_number = 1.
dif_phase_v_h_deg = 0.

""" Animation control """
is_play = False

""" Axis vectors """
vector_x_axis = np.array([1., 0., 0.])
vector_y_axis = np.array([0., 1., 0.])
vector_z_axis = np.array([0., 0., 1.])

""" Create figure and axes """
title_tk = "Möbius strip"
title_ax0 = title_tk

x_min = -2.
x_max = 2.
y_min = -2.
y_max = 2.
z_min = -2.
z_max = 2.

fig = Figure(facecolor='black')
ax0 = fig.add_subplot(111, projection="3d")
ax0.set_box_aspect((4, 4, 4))
ax0.grid()
ax0.set_title(title_ax0, color="white")
ax0.set_xlabel("x")
ax0.set_ylabel("y")
ax0.set_zlabel("z")
ax0.set_xlim(x_min, x_max)
ax0.set_ylim(y_min, y_max)
ax0.set_zlim(z_min, z_max)

ax0.set_facecolor("black")
ax0.axis('off')

x_min = -2.
x_max = 2.
y_min = -2.
y_max = 2.
z_min = -2.
z_max = 2.


""" Embed in Tkinter """
root = tk.Tk()
root.title(title_tk)
canvas = FigureCanvasTkAgg(fig, root)
canvas.get_tk_widget().pack(expand=True, fill="both")

toolbar = NavigationToolbar2Tk(canvas, root)
canvas.get_tk_widget().pack()

""" Global objects of Tkinter """
var_num_lines = tk.StringVar(root)
var_is_spiral = tk.IntVar(root)
var_pitch_spiral = tk.StringVar(root)
var_num_lap = tk.StringVar(root)
var_num_twist = tk.StringVar(root)
var_is_wave = tk.IntVar(root)
var_amp_v = tk.StringVar(root)
var_amp_h = tk.StringVar(root)
var_wn = tk.StringVar(root)
var_dif_phase_deg = tk.StringVar(root)

""" Classes and functions """


class Counter:
    def __init__(self, is3d=None, ax=None, xy=None, z=None, label=""):
        self.is3d = is3d if is3d is not None else False
        self.ax = ax
        self.x, self.y = xy[0], xy[1]
        self.z = z if z is not None else 0
        self.label = label

        self.count = 0

        if not is3d:
            self.txt_step = self.ax.text(self.x, self.y, self.label + str(self.count), color="white")
        else:
            self.txt_step = self.ax.text2D(self.x, self.y, self.label + str(self.count), color="white")
            self.xz, self.yz, _ = proj3d.proj_transform(self.x, self.y, self.z, self.ax.get_proj())
            self.txt_step.set_position((self.xz, self.yz))

    def count_up(self):
        self.count += 1
        self.txt_step.set_text(self.label + str(self.count))

    def reset(self):
        self.count = 0
        self.txt_step.set_text(self.label + str(self.count))

    def get(self):
        return self.count


class CrossLine:
    def __init__(self, ax, size, radius, phi, amplitude_v, amplitude_h):
        self.ax = ax
        self.size = size
        self.radius = radius
        self.phi = phi
        self.amplitude_v = amplitude_v
        self.amplitude_h = amplitude_h

        self.pitch = 0.

        self.is_wave = False
        self.wave_number = 1.
        self.dif_phase_v_h = 0.
        self.progress_phase = 0.

        self.origin = np.array([np.cos(self.phi), np.sin(self.phi), 0.]) * self.radius
        self.vertical_axis = np.array([0., 0., 1.])
        self.horizontal_axis = np.array([1., 0., 0.])

        rot_matrix_z = Rotation.from_rotvec(self.phi * vector_z_axis)
        self.horizontal_axis_rotated = rot_matrix_z.apply(self.horizontal_axis)
        self.vertical_axis_rotated = self.vertical_axis

        line_v0 = zip(self.origin, self.vertical_axis_rotated * self.amplitude_v + self.origin)
        self.plt_line_v0, = self.ax.plot(*line_v0, linewidth=1, linestyle="-", color="red")

        line_v1 = zip(self.origin, - self.vertical_axis_rotated * self.amplitude_v + self.origin)
        self.plt_line_v1, = self.ax.plot(*line_v1, linewidth=1, linestyle="-", color="green")

        line_h0 = zip(self.origin, self.horizontal_axis_rotated * self.amplitude_h + self.origin)
        self.plt_line_h0, = self.ax.plot(*line_h0, linewidth=0.5, linestyle=":", color="blue")

        line_h1 = zip(self.origin, - self.horizontal_axis_rotated * self.amplitude_h + self.origin)
        self.plt_line_h1, = self.ax.plot(*line_h1, linewidth=0.5, linestyle=":", color="darkorange")

    def remove(self):
        self.plt_line_v0.remove()
        self.plt_line_v1.remove()
        self.plt_line_h0.remove()
        self.plt_line_h1.remove()

    def set_radius(self, radius):
        self.radius = radius
        self.update_diagrams()

    def set_pitch(self, pitch):
        self.pitch = pitch
        self.update_diagrams()

    def update_diagrams(self):
        self.origin = (np.array([np.cos(self.phi), np.sin(self.phi), 0.]) *
                       (self.radius + self.pitch * self.phi / (2. * np.pi)))

        if self.is_wave:
            if np.cos(self.phi * self.wave_number + self.progress_phase) >= 0:
                amp_v0 = np.cos(self.phi * self.wave_number + self.progress_phase) * self.amplitude_v
                amp_v1 = 0.
            else:
                amp_v0 = 0.
                amp_v1 = np.abs(np.cos(self.phi * self.wave_number + self.progress_phase) * self.amplitude_v)

            if np.cos(self.phi * self.wave_number + self.progress_phase) >= 0:
                amp_h0 = np.cos((self.phi + self.dif_phase_v_h) * self.wave_number + self.progress_phase) * self.amplitude_h
                amp_h1 = 0.
            else:
                amp_h0 = 0.
                amp_h1 = np.abs(np.cos((self.phi + self.dif_phase_v_h) * self.wave_number + self.progress_phase) * self.amplitude_h)

        else:
            amp_v0 = self.amplitude_v
            amp_v1 = self.amplitude_v
            amp_h0 = self.amplitude_h
            amp_h1 = self.amplitude_h

        line_v0 = zip(self.origin, self.vertical_axis_rotated * amp_v0 + self.origin)
        self.plt_line_v0.set_data_3d(*line_v0)

        line_v1 = zip(self.origin, - self.vertical_axis_rotated * amp_v1 + self.origin)
        self.plt_line_v1.set_data_3d(*line_v1)

        line_h0 = zip(self.origin, self.horizontal_axis_rotated * amp_h0 + self.origin)
        self.plt_line_h0.set_data_3d(*line_h0)

        line_h1 = zip(self.origin, - self.horizontal_axis_rotated * amp_h1 + self.origin)
        self.plt_line_h1.set_data_3d(*line_h1)

    def roll(self, angle):
        self.vertical_axis = np.array([0., 0., 1.])
        self.horizontal_axis = np.array([1., 0., 0.])
        self.vertical_axis_rotated = self.vertical_axis

        rot_matrix_z = Rotation.from_rotvec(self.phi * vector_z_axis)
        self.horizontal_axis_rotated = rot_matrix_z.apply(self.horizontal_axis)

        center_axis = np.cross(self.vertical_axis_rotated, self.horizontal_axis_rotated)
        rot_matrix_center = Rotation.from_rotvec(angle * center_axis)
        self.vertical_axis_rotated = rot_matrix_center.apply(self.vertical_axis_rotated)
        self.horizontal_axis_rotated = rot_matrix_center.apply(self.horizontal_axis_rotated)

        self.update_diagrams()

    def set_amplitude_v(self, value):
        self.amplitude_v = value
        self.update_diagrams()

    def set_amplitude_h(self, value):
        self.amplitude_h = value
        self.update_diagrams()

    def set_is_wave(self, value):
        self.is_wave = value
        self.update_diagrams()

    def set_wave_number(self, value):
        self.wave_number = value
        self.update_diagrams()

    def set_dif_phase_deg(self, value):
        self.dif_phase_v_h = np.deg2rad(value)
        self.update_diagrams()

    def set_pregress_phase(self, value):
        self.progress_phase = value
        self.update_diagrams()


class Path:
    def __init__(self, ax, line_width, color):
        self.ax = ax
        self.line_width = line_width
        self.color = color

        self.is_draw_path = False

        self.x_path = []
        self.y_path = []
        self.z_path = []
        self.path, = self.ax.plot(np.array(self.x_path), np.array(self.y_path), np.array(self.z_path),
                                  color=self.color, linewidth=self.line_width)

    def append_path(self, position):
        if self.is_draw_path:
            self.x_path.append(position[0])
            self.y_path.append(position[1])
            self.z_path.append(position[2])
            self.update_path()

    def update_path(self):
        self.path.set_data_3d(np.array(self.x_path), np.array(self.y_path), np.array(self.z_path))

    def clear_path(self):
        self.x_path = []
        self.y_path = []
        self.z_path = []
        self.update_path()

    def set_is_draw_path(self, value):
        self.is_draw_path = value


def update_diagram():
    phase = np.deg2rad(cnt.get())
    for line in cross_lines:
        line.set_pregress_phase(phase)


def set_is_spiral_base(value):
    if value:
        for i in range(len(cross_lines)):
            cross_lines[i].set_pitch(pitch_spiral)
    else:
        for i in range(len(cross_lines)):
            cross_lines[i].set_pitch(0.)


def set_cross_lines(num):
    global num_cross_lines, cross_lines
    num_cross_lines = num

    for line in cross_lines:
        line.remove()
    cross_lines.clear()

    for i in range(num_cross_lines):
        phi = i * (2. * np.pi * num_lap) / num_cross_lines
        cross_line = CrossLine(ax0, 1, 1, phi, amplitude_v, amplitude_h)
        cross_line.roll(phi * num_twist)
        cross_lines.append(cross_line)

    set_pitch_spiral(pitch_spiral)
    set_is_wave(var_is_wave.get())


def set_pitch_spiral(value):
    global pitch_spiral
    pitch_spiral = value

    for i in range(num_cross_lines):
        if var_is_spiral.get():
            cross_lines[i].set_pitch(pitch_spiral)
        else:
            cross_lines[i].set_pitch(0.)


def set_num_lap(value):
    global num_lap
    num_lap = value

    set_cross_lines(num_cross_lines)


def set_num_twist(value):
    global num_twist
    num_twist = value

    for i in range(num_cross_lines):
        phi = i * (2. * np.pi * num_lap) / num_cross_lines
        cross_lines[i].roll(phi * num_twist)


def set_is_wave(value):
    for line in cross_lines:
        line.set_is_wave(value)


def set_amplitude_v(value):
    for line in cross_lines:
        line.set_amplitude_v(value)


def set_amplitude_h(value):
    for line in cross_lines:
        line.set_amplitude_h(value)


def set_wave_number(value):
    for line in cross_lines:
        line.set_wave_number(value)


def set_dif_phase_deg(value):
    for line in cross_lines:
        line.set_dif_phase_deg(value)


def create_parameter_setter():
    # Amplitude
    frm_amp = ttk.Labelframe(root, relief="ridge", text="Amplitude", labelanchor="n")
    frm_amp.pack(side='left', fill=tk.Y)

    lbl_amp_v = tk.Label(frm_amp, text="Vertical")
    lbl_amp_v.pack(side='left')

    # var_amp_v = tk.StringVar(root)
    var_amp_v.set(str(amplitude_v))
    spn_amp_v = tk.Spinbox(
        frm_amp, textvariable=var_amp_v, format="%.1f", from_=0., to=5., increment=0.1,
        command=lambda: set_amplitude_v(float(var_amp_v.get())), width=5
    )
    spn_amp_v.pack(side="left")

    lbl_amp_h = tk.Label(frm_amp, text="Horizontal")
    lbl_amp_h.pack(side='left')

    # var_amp_h = tk.StringVar(root)
    var_amp_h.set(str(amplitude_v))
    spn_amp_h = tk.Spinbox(
        frm_amp, textvariable=var_amp_h, format="%.1f", from_=0., to=5., increment=0.1,
        command=lambda: set_amplitude_h(float(var_amp_h.get())), width=5
    )
    spn_amp_h.pack(side="left")

    # Number of lines
    frm_num = ttk.Labelframe(root, relief="ridge", text="Number of lines", labelanchor='n')
    frm_num.pack(side="left", fill=tk.Y)

    # var_num_lines = tk.StringVar(root)
    var_num_lines.set(str(num_cross_lines))
    spn_num_lines = tk.Spinbox(
        frm_num, textvariable=var_num_lines, format="%.0f", from_=1, to=360, increment=1,
        command=lambda: set_cross_lines(int(var_num_lines.get())), width=5
    )
    spn_num_lines.pack(side="left")

    # Spiral of base
    frm_spiral = ttk.Labelframe(root, relief="ridge", text="Spiral", labelanchor="n")
    frm_spiral.pack(side='left', fill=tk.Y)

    # var_is_spiral = tk.IntVar(root)
    chk_is_spiral = tk.Checkbutton(frm_spiral, text="Apply", variable=var_is_spiral,
                                   command=lambda: set_is_spiral_base(var_is_spiral.get()))
    chk_is_spiral.pack(side='left')
    var_is_spiral.set(False)

    # var_pitch_spiral = tk.StringVar(root)
    var_pitch_spiral.set(str(pitch_spiral))
    spn_pitch_spiral = tk.Spinbox(
        frm_spiral, textvariable=var_pitch_spiral, format="%.1f", from_=0.0, to=2.0, increment=0.1,
        command=lambda: set_pitch_spiral(float(var_pitch_spiral.get())), width=5
    )
    spn_pitch_spiral.pack(side="left")

    # Number of lap
    frm_lap = ttk.Labelframe(root, relief="ridge", text="Number of lap", labelanchor="n")
    frm_lap.pack(side='left', fill=tk.Y)

    # var_num_lap = tk.StringVar(root)
    var_num_lap.set(str(num_lap))
    spn_num_lap = tk.Spinbox(
        frm_lap, textvariable=var_num_lap, format="%.1f", from_=0.1, to=6.0, increment=0.1,
        command=lambda: set_num_lap(float(var_num_lap.get())), width=5
    )
    spn_num_lap.pack(side="left")

    # Number of twist
    frm_twist = ttk.Labelframe(root, relief="ridge", text="Number of twist", labelanchor="n")
    frm_twist.pack(side='left', fill=tk.Y)

    # var_num_twist = tk.StringVar(root)
    var_num_twist.set(str(num_twist))
    spn_num_twist = tk.Spinbox(
        frm_twist, textvariable=var_num_twist, format="%.1f", from_=0., to=10., increment=0.5,
        command=lambda: set_num_twist(float(var_num_twist.get())), width=5
    )
    spn_num_twist.pack(side="left")

    # Wave
    frm_wave = ttk.Labelframe(root, relief="ridge", text="Wave", labelanchor="n")
    frm_wave.pack(side='left', fill=tk.Y)

    # var_is_wave = tk.IntVar(root)
    chk_is_wave = tk.Checkbutton(frm_wave, text="Apply", variable=var_is_wave,
                                 command=lambda: set_is_wave(var_is_wave.get()))
    chk_is_wave.pack(side='left')
    var_is_wave.set(False)

    lbl_wave_num = tk.Label(frm_wave, text="Wave number")
    lbl_wave_num.pack(side='left')

    # var_wn = tk.StringVar(root)
    var_wn.set(str(wave_number))
    spn_wn = tk.Spinbox(
        frm_wave, textvariable=var_wn, format="%.0f", from_=1, to=10, increment=1,
        command=lambda: set_wave_number(float(var_wn.get())), width=5
    )
    spn_wn.pack(side="left")

    lbl_dif = tk.Label(frm_wave, text="Diff. phase between V, H")
    lbl_dif.pack(side='left')

    # var_dif_phase_deg = tk.StringVar(root)
    var_dif_phase_deg.set(str(dif_phase_v_h_deg))
    spn_dif_phase_deg = tk.Spinbox(
        frm_wave, textvariable=var_dif_phase_deg, format="%.0f", from_=0, to=360, increment=1,
        command=lambda: set_dif_phase_deg(float(var_dif_phase_deg.get())), width=5
    )
    spn_dif_phase_deg.pack(side="left")


def create_animation_control():
    frm_anim = ttk.Labelframe(root, relief="ridge", text="Animation", labelanchor="n")
    frm_anim.pack(side="left", fill=tk.Y)
    btn_play = tk.Button(frm_anim, text="Play/Pause", command=switch)
    btn_play.pack(side="left")
    btn_reset = tk.Button(frm_anim, text="Reset", command=reset)
    btn_reset.pack(side="left")
    # btn_clear = tk.Button(frm_anim, text="Clear path", command=lambda: aaa())
    # btn_clear.pack(side="left")


def create_center_lines():
    ln_axis_x = art3d.Line3D([x_min, x_max], [0., 0.], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax0.add_line(ln_axis_x)
    ln_axis_y = art3d.Line3D([0., 0.], [y_min, y_max], [0., 0.], color="gray", ls="-.", linewidth=1)
    ax0.add_line(ln_axis_y)
    ln_axis_z = art3d.Line3D([0., 0.], [0., 0.], [z_min, z_max], color="gray", ls="-.", linewidth=1)
    ax0.add_line(ln_axis_z)


def draw_static_diagrams():
    create_center_lines()


def reset():
    global is_play
    if is_play:
        is_play = not is_play
    cnt.reset()
    update_diagram()


def switch():
    global is_play
    if is_play:
        is_play = False
    else:
        is_play = True


def update(f):
    if is_play:
        cnt.count_up()
        update_diagram()


""" main loop """
if __name__ == "__main__":
    cnt = Counter(ax=ax0, is3d=True, xy=np.array([x_min, y_max]), z=z_max, label="Step=")
    draw_static_diagrams()
    create_animation_control()
    create_parameter_setter()

    for i_ in range(num_cross_lines):
        cross_line = CrossLine(ax0, 1, 1, i_ * (2. * np.pi * num_lap) / num_cross_lines,
                               amplitude_v, amplitude_h)
        cross_lines.append(cross_line)

    set_num_twist(num_twist)

    # ax0.legend(loc='lower right', fontsize=8)

    anim = animation.FuncAnimation(fig, update, interval=100, save_count=100)
    root.mainloop()
