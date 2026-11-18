import logging

logger = logging.getLogger(__name__)


def check_plugins(config_parser, head):
    """Check if the plugins exist in vcf file

    Args:
        config_parser (ConfigObj): A config object with the plugins
        head (HeaderParser): A vcf header object

    Returns:
        bool: If all tests passed or not
    """
    all_pass = True
    for plugin in config_parser.plugins:
        plugin_object = config_parser.plugins[plugin]
        logger.debug("Checking plugin {0}".format(plugin))
        if plugin_object.field == "INFO":
            info_key = plugin_object.info_key
            if info_key not in head.info_dict:
                logger.warning(
                    "INFO field {0} is not in vcf INFO. This field will not be scored.".format(
                        info_key
                    )
                )
                all_pass = False
            else:
                logger.debug("INFO field {0} was found in vcf INFO.".format(info_key))
            if info_key in ("CSQ", "ANN"):
                csq_key = plugin_object.csq_key
                if csq_key not in head.annotation_columns(info_key):
                    logger.warning(
                        "%s field %s is not in %s annotation. This field will not be scored.",
                        info_key,
                        csq_key,
                        info_key,
                    )
                    all_pass = False
                else:
                    logger.debug(
                        "%s field %s was found in %s annotation.", info_key, csq_key, info_key
                    )
    return all_pass
