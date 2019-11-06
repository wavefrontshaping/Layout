import logging


def get_logger(name):
        loglevel = logging.DEBUG
        logger = logging.getLogger(name)
        if not getattr(logger, 'handler_set', None):
            logger.setLevel(logging.INFO)
            logFormatter = logging.Formatter("%(asctime)s - %(name)-10.10s [%(levelname)-7.7s]  %(message)s") #[%(threadName)-12.12s] 
            consoleHandler = logging.StreamHandler()
            consoleHandler.setFormatter(logFormatter)
            logger.addHandler(consoleHandler)
            logger.setLevel(loglevel)
            logger.handler_set = True
        return logger
    