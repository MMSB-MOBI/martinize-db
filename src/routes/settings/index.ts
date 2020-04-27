import { Router } from 'express';
import { promises as FsPromise } from 'fs';
import { SETTINGS_FILE } from '../../constants';
import logger from '../../logger';
import { errorCatcher, methodNotAllowed } from '../../helpers';
import Errors, { ErrorType } from '../../Errors';
import { SettingsJson } from '../../types';
import { Database } from '../../Entities/CouchHelper';
import MembraneBuilder from '../../Builders/MembraneBuilder';
import RadiusDatabase from '../../Entities/RadiusDatabase';

const SettingsRouter = Router();

// Settings file "settings.json" at project root
SettingsRouter.get('/', (_, res) => {
  FsPromise.readFile(SETTINGS_FILE, "utf-8")
    .then(file => {
      res.json(JSON.parse(file));
    })
    .catch(e => {
      logger.error("Unable to get settings file.", e);
      Errors.throw(ErrorType.Server);
    })
    .catch(errorCatcher(res));
});

SettingsRouter.post('/', (req, res) => {
  if (req.full_user?.role !== "admin") {
    Errors.throw(ErrorType.Forbidden);
  } 

  if (typeof req.body !== 'object') {
    return Errors.throw(ErrorType.Format);
  }

  (async () => {
    const file: SettingsJson = JSON.parse(await FsPromise.readFile(SETTINGS_FILE, 'utf-8'));
    const new_file = { ...file };

    for (const prop in req.body) {
      // @ts-ignore
      new_file[prop] = req.body[prop];
    }

    await FsPromise.writeFile(SETTINGS_FILE, JSON.stringify(new_file, null, 2));

    res.json(new_file);
  })().catch(errorCatcher(res));
});

SettingsRouter.all('/', methodNotAllowed(['GET', 'POST']));

SettingsRouter.get('/lipids', async (req, res) => {
  let ff: string | undefined = undefined;

  if (req.query.force_field) {
    ff = req.query.force_field as string;
  }

  /// WITH DATABASE
  // res.json(await Database.lipid.getAvailableNames(ff));
  /// WITH FILES
  if (ff) {
    const prefix = RadiusDatabase.FORCE_FIELD_TO_MARTINI_VERSION[ff];

    if (!prefix) {
      return [];
    }

    res.json(MembraneBuilder.SUPPORTED_LIPIDS[prefix] ?? []);
  }
  else {
    const lipids = new Set<string>();

    for (const prefix in MembraneBuilder.SUPPORTED_LIPIDS) {
      for (const lipid of MembraneBuilder.SUPPORTED_LIPIDS[prefix]) {
        lipids.add(lipid);
      }
    }

    res.json([...lipids]);
  }
});

SettingsRouter.all('/lipids', methodNotAllowed(['GET']));

export default SettingsRouter;
