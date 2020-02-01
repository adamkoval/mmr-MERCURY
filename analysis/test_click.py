import matplotlib.pyplot as plt
import numpy as np

def onpick(event):
    point = event.artist
#    global x, y
    x, y = point.get_data()

    global coords
    coords.append((x, y))

    print('x, y = {}, {}'.format(x, y))

if __name__ == "__main__":
    fig = plt.figure()
    for i in range(100):
        randx = np.random.uniform(0, 1)
        randy = np.random.uniform(0, 1)
        plt.plot(randx, randy, 'ko', picker=4)
#    fig.canvas.mpl_connect('key_press_event', ontype)
#    fig.canvas.mpl_connect('button_press_event', onclick)
#    fig.canvas.mpl_connect('button_release_event', onclick)
    coords = []
    fig.canvas.mpl_connect('pick_event', onpick)
    plt.show()
