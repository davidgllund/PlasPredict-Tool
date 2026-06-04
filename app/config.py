"""
Configuration for Plasmid Host Range Predictor Web Application
"""
import os
from pathlib import Path

# Project root directory
PROJECT_ROOT = Path(__file__).parent


class Config:
    """Base configuration"""
    
    # Flask settings
    FLASK_ENV = 'development'
    DEBUG = False
    TESTING = False
    
    # Security
    SECRET_KEY = os.environ.get('SECRET_KEY', 'dev-secret-key-change-in-production')
    SESSION_COOKIE_SECURE = True
    SESSION_COOKIE_HTTPONLY = True
    SESSION_COOKIE_SAMESITE = 'Lax'
    
    # File upload
    MAX_CONTENT_LENGTH = 50 * 1024 * 1024  # 50MB
    UPLOAD_FOLDER = os.path.join(PROJECT_ROOT, 'uploads')
    ALLOWED_EXTENSIONS = {'fna', 'fasta', 'fa', 'txt'}
    
    # Model paths - support both development and container environments
    MODEL_PATH = os.environ.get(
        'MODEL_PATH',
        PROJECT_ROOT / 'models' / 'plaspredict_model.pkl'
    )
    HMM_PATH = os.environ.get(
        'HMM_PATH',
        PROJECT_ROOT / 'models' / 'conjscan_models'
    )
    
    PLASMIDFINDER_DB_PATH = os.environ.get(
        'PLASMIDFINDER_DB_PATH',
        os.path.expanduser('$HOME/app/db/plasmidfinder/database/')
    )
    
    # Logging
    LOG_LEVEL = os.environ.get('LOG_LEVEL', 'DEBUG')  # Use DEBUG in Docker
    LOG_FILE = os.path.join(PROJECT_ROOT, 'logs', 'app.log')


class DevelopmentConfig(Config):
    """Development configuration"""
    FLASK_ENV = 'development'
    DEBUG = True
    SESSION_COOKIE_SECURE = False


class ProductionConfig(Config):
    """Production configuration"""
    FLASK_ENV = 'production'
    DEBUG = False
    TESTING = False
    
    # Require secure cookies in production (behind proxy with HTTPS)
    SESSION_COOKIE_SECURE = False  # Set to False when behind reverse proxy
    
    # Stricter security
    PREFERRED_URL_SCHEME = 'https'


class TestingConfig(Config):
    """Testing configuration"""
    TESTING = True
    DEBUG = True
    
    # Use in-memory temp directory for testing
    UPLOAD_FOLDER = '/tmp/plasmid-predictor-test'
    
    # Disable CSRF for testing
    WTF_CSRF_ENABLED = False


# Load configuration based on environment
config = {
    'development': DevelopmentConfig,
    'production': ProductionConfig,
    'testing': TestingConfig,
    'default': DevelopmentConfig
}

def get_config(env=None):
    """Get configuration object"""
    if env is None:
        env = os.environ.get('FLASK_ENV', 'development')
    return config.get(env, config['default'])
