"""
## Schema


THE FOLLOWING ARE _NOT_ IN ANY CONFIG.YAML FILE!!!
it would be confusing to the user to have multiple sources of truth


fireworks_config_mode: user (.fireworks) | env | manual (site-packages/fireworks)
fireworks_config_dir: (you pick), defaults to whatever fireworks_config_mode is set to

It makes, if requested:
(Fireworks)
- FW_config.yaml
- my_launchpad.yaml
    - Call fireworks.scripts.lpad_run.init_yaml?
    - or just ask for database string
- my_fworker.yaml
    - has to point to db.json
    - needs the vasp CMD
- my_qadapter.yaml
    - activate_env_script
(Atomate)
- db.json
    - This could be anywhere at all, atomate only knows from my_fworker
(PyMatGen)
- .pmgrc.yaml
"""

import os
import json

import click
import yaml


INIT_DEFAULT_SETTINGS = {

}


def cmd_easy_init():

    # https://eng.localytics.com/exploring-cli-best-practices/
    # https://click.palletsprojects.com/en/7.x/options/#option-prompting
    # https://click.palletsprojects.com/en/7.x/prompts/#option-prompts
    # https://click.palletsprojects.com/en/7.x/quickstart/#screencast-and-examples
    # https://click.palletsprojects.com/en/7.x/setuptools/#setuptools-integration

    # TODO: This is temporarily hardcoded
    config = dict()

    config['fireworks_config_dir'] = os.path.expanduser('~/.fireworks')
    config['atomate_config_dir'] = os.path.expanduser('~/.atomate')

    config['db_json_path'] = os.path.join(config['atomate_config_dir'], 'db.json')
    config['db_json_config'] = {
    }

    config['launchpad_path'] = os.path.join(config['fireworks_config_dir'], 'my_launchpad.yaml')
    config['launchpad_config'] = {
        'host': 'mongodb://fw_BL_931-2_CRN_admin:ww11122wfg64bwww@mongodb07.nersc.gov:27017/fw_BL_931-2_CRN',
    }

    config['fworker_path'] = os.path.join(config['fireworks_config_dir'], 'my_fworker.yaml')
    config['fworker_config'] = {
        'db_file': config['db_json_path'],
    }

    config['qadapter_path'] = os.path.join(config['fireworks_config_dir'], 'my_qadapter.yaml')
    config['qadapter_config'] = {
    }

    write_out(config)


def write_out(config):

    # Make the proper paths if they don't already exist
    # TODO: Make sure that I can write to launchpad_path and whatnot
    if not os.path.exists(config['fireworks_config_dir']):
        os.mkdir(config['fireworks_config_dir'])
    if not os.path.exists(config['atomate_config_dir']):
        os.mkdir(config['atomate_config_dir'])

    # my_launchpad.yaml
    with open(config['launchpad_path'], 'w') as file:
        file.write(yaml.dump(config['launchpad_config']))

    # my_fworker.yaml
    with open(config['fworker_path'], 'w') as file:
        file.write(yaml.dump(config['fworker_config']))

    # my_qadapter.yaml
    with open(config['qadapter_path'], 'w') as file:
        file.write(yaml.dump(config['qadapter_config']))

    # db.json
    json.dump(config['db_json_config'], config['db_json_path'])

    # .pmgrc.yaml