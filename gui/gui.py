import tkinter as tk

# #                     # #
#   General use classes   #
# #                     # #
class GenericPage(tk.Frame):
    def __init__(self, parent, controller, page_name):
        all_buttons = {"Home": 'HomePage',
                "Setup": 'SetupPage',
                "Simulation": 'SimPage',
                "Analysis": 'AnalysisPage'}
        buttons = {key: all_buttons[key] for key in all_buttons if key != page_name}
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


class ParamEditor(tk.Toplevel):
    def __init__(self, parent, file):
        tk.Toplevel.__init__(self, parent)
        text_window = TextWindow(self, file)
        save_button = tk.Button(self, text="Save changes",
                command=lambda: text_window.save_file())
        save_button.pack()


class TextWindow(tk.Frame):
    def __init__(self, parent, file):
        tk.Frame.__init__(self, parent)
        self.file = file
        self.parent = parent

        text = open(self.file, 'r').read()
        self.textbox = tk.Text(self.parent)
        self.textbox.pack()
        self.textbox.insert(1.0, text)

    def save_file(self):
        text = self.textbox.get(1.0, tk.END)
        f = open(self.file, 'w')
        f.write(text)
        f.close()
        self.parent.destroy()


# #       # #
#   Pages   #
# #       # #
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
        page = self.pages[page_name]
        page.tkraise()


class HomePage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Home")


class SetupPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Setup")
        paramin = "text.txt"
        paramin_button = tk.Button(self, text="Edit param.in",
                command=lambda: Popup(self, paramin))
        paramin_button.pack()


class AnalysisPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Analysis")


class SimPage(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        GenericPage(self, controller, "Simulation")


if __name__=="__main__":
    app = MainApp()
    app.mainloop()
