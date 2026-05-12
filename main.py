import tkinter as tk
from tkinter import ttk, scrolledtext
import traceback
import sys

def show_copyable_error(title, error_msg):
    """에러 내용을 복사할 수 있는 창을 띄웁니다."""
    root = tk.Tk()
    root.title(title)
    root.geometry("600x400")
    
    ttk.Label(root, text="An error occurred. You can copy the details below:", padding=10).pack(anchor="w")
    
    # 텍스트 복사가 가능한 창
    txt = scrolledtext.ScrolledText(root, padx=10, pady=10)
    txt.pack(expand=True, fill="both", padx=10, pady=5)
    txt.insert(tk.END, error_msg)
    txt.config(state="disabled") # 읽기 전용이지만 복사는 가능
    
    btn_f = ttk.Frame(root)
    btn_f.pack(fill="x", pady=10)
    
    def copy_to_clipboard():
        root.clipboard_clear()
        root.clipboard_append(error_msg)
        root.update()
        
    ttk.Button(btn_f, text="Copy to Clipboard", command=copy_to_clipboard).pack(side="left", padx=20)
    ttk.Button(btn_f, text="Close", command=root.destroy).pack(side="right", padx=20)
    
    root.mainloop()

try:
    from app import YeastApp
except Exception:
    show_copyable_error("Startup Error", traceback.format_exc())
    sys.exit(1)

if __name__ == "__main__":
    try:
        root = tk.Tk()
        app = YeastApp(root)
        root.mainloop()
    except Exception:
        show_copyable_error("Runtime Error", traceback.format_exc())
