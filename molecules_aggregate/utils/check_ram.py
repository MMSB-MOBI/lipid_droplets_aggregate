from multiprocessing import Process, Manager
import psutil, logging, time

def get_ram():
    while True:
        ram = psutil.virtual_memory()._asdict()
        used = ram["used"]
        available = ram["available"]
        if available - used < 10 * 1024 * 1024 * 1024 : 
            logging.warning("Less than 10 Go RAM available !")
            if available - used < 5 * 1024 * 1024 * 1024 : 
                logging.error("Too much RAM used. Less than 5 Go available")
                exit()
        if used > MAX_RAM.value:
            MAX_RAM.value = used
        time.sleep(1)

def ram_human_readable(ram):
    units = ["o", "Ko", "Mo", "Go", "To"]
    mem = float(ram)
    unit = 0
    while mem >= 1024:
        mem = mem / 1024
        unit += 1
    return f"{round(mem,3)} {units[unit]}"

def start():
    global RAM_PROCESS
    global MAX_RAM

    MAX_RAM = Manager().Value("i", 0)
    RAM_PROCESS =  Process(target = get_ram)
    RAM_PROCESS.start()

def stop():
    RAM_PROCESS.terminate()
    logging.info(f"MAX RAM USED : {ram_human_readable(MAX_RAM.value)}")