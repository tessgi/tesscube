import configparser
import os

from . import _package_cache_dir, log

CONFIG_FILE_PATH = os.path.join(_package_cache_dir, "config.ini")


def make_default_config():
    config = configparser.ConfigParser()
    config["Settings"] = {}
    log.debug("Writing default config file.")
    with open(CONFIG_FILE_PATH, "w") as f:
        config.write(f)
    log.debug("Writen.")


def load_config():
    config = configparser.ConfigParser()
    if not os.path.exists(CONFIG_FILE_PATH):
        log.debug("No configuration file found, creating the default file.")
        make_default_config()
    else:
        config.read(CONFIG_FILE_PATH)
    return config


def save_config(config):
    with open(CONFIG_FILE_PATH, "w") as f:
        config.write(f)


def set_config_parameter(parameter_name: str, parameter_value: str):
    """Update a configuration parameter"""
    config = load_config()
    config["Settings"][parameter_name] = f"{parameter_value}"
    save_config(config)


def get_config_parameter(parameter_name: str):
    return load_config()["Settings"][parameter_name]
