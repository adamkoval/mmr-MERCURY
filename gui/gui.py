import tkinter as tk

class MainApp(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)

        mainframe = tk.Frame(self)
        mainframe.pack(side="top", fill="both", expand=True)
        mainframe.grid_rowconfigure(0, weight=1)
        mainframe.grid_columnconfigure(0, weight=1)

        self.pages = {}
        for _page in (HomePage, SimPage, AnalysisPage, SetupPage):
            page_name = _page.__name__
            page = _page(parent=mainframe, controller=self)
            self.pages[page_name] = page
            page.grid(row=0, column=0, sticky="NESW")

        self.show_page("HomePage")

    def show_page(self, page_name):
        print("Raising page {}".format(page_name)) # DEBUG
        page = self.pages[page_name]
        page.tkraise()


class SetupPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Setup")


class AnalysisPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Analysis")


class SimPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Simulation")


class HomePage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Home")


class GenericPage(tk.Frame):
    def __init__(self, parent, controller, page_name):
        all_buttons = {"Home": 'HomePage',
                "Setup": 'SetupPage',
                "Simulation": 'SimPage',
                "Analysis": 'AnalysisPage'}
        buttons = {button: all_buttons[button] for button in all_buttons if button != page_name}
        tk.Frame.__init__(self, parent)
        label = tk.Label(parent, text=page_name)
        label.pack()
        navbar = NavBar(parent, controller, buttons)


class NavButton(tk.Frame):
    def __init__(self, parent, controller, text, dest, col):
        tk.Frame.__init__(self, parent)
        button = tk.Button(parent, text=text,
                command=lambda: controller.show_page(dest))
        button.grid(row=0, column=col)


class NavBar(tk.Frame):
    def __init__(self, parent, controller, buttons):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        for i, text in enumerate(buttons):
            dest = buttons[text]
            NavButton(self, controller, text, dest, i)
        self.pack()


if __name__=="__main__":
    app = MainApp()
    app.mainloop()
